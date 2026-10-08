#!/usr/bin/env python3
"""Selective CI: change detection and per-area fingerprints.

Usage in CI::

    python3 ci/fingerprint.py --base "$BASE" --head HEAD >> "$GITHUB_OUTPUT"
    python3 ci/fingerprint.py --decide >> "$GITHUB_OUTPUT"
    python3 ci/fingerprint.py --check-workflow .github/workflows/ci.yml

Every doubt resolves to "run the job": parse failures fall back to raw bytes,
git failures mark every area affected, and an unusable fingerprint gets a
random key so it can never match a recorded pass marker.

Standard library only; keep compatible with Python 3.9.
"""

import argparse
import hashlib
import io
import json
import os
import re
import subprocess
import sys
import tokenize
import uuid

SCRIPT_VERSION = "selective-ci-1"
HERE = os.path.dirname(os.path.abspath(__file__))
CONFIG = os.path.join(HERE, "areas.json")
WORKFLOW = ".github/workflows/ci.yml"
EMPTY_TREE = "4b825dc642cb6eb9a060e54bf8d69288fbee4904"


class Fallback(Exception):
    """Raised by a normalizer that cannot be sure; caller uses raw bytes."""


# --------------------------------------------------------------------------
# Globs

def glob_to_regex(pattern):
    """Translate a path glob to a regex: ** spans directories, * and ? do not."""
    out, i, n = [], 0, len(pattern)
    while i < n:
        c = pattern[i]
        if pattern.startswith("**/", i):
            out.append("(?:.*/)?")
            i += 3
        elif pattern.startswith("**", i):
            out.append(".*")
            i += 2
        elif c == "*":
            out.append("[^/]*")
            i += 1
        elif c == "?":
            out.append("[^/]")
            i += 1
        else:
            out.append(re.escape(c))
            i += 1
    return re.compile("".join(out) + r"\Z")


class Matcher:
    def __init__(self, patterns):
        self.regexes = [glob_to_regex(p) for p in patterns]

    def __call__(self, path):
        return any(r.match(path) for r in self.regexes)


# --------------------------------------------------------------------------
# Normalization

DIRECTIVES = re.compile(
    r"^#!|-\*-\s*coding|coding[:=]|#\s*type:|noqa|pragma|pylint:|mypy:|pyright:|"
    r"\bfmt:|isort:|yamllint|shellcheck|"
    r"^//go:|^//\s*\+build|@ts-|^///\s*<reference|eslint-|prettier-ignore|"
    r"istanbul ignore|c8 ignore|[#@]__PURE__|@license|@preserve|#(?:end)?region|"
    r"NOLINT|clang-format|LCOV_EXCL|cppcheck-suppress"
)

C_LIKE = {".c", ".h", ".C", ".H", ".cpp", ".hpp", ".cc", ".hh", ".cxx", ".hxx",
          ".ipp", ".tpp", ".inl"}
JS_LIKE = {".js", ".jsx", ".ts", ".tsx", ".mjs", ".cjs"}
TRIPLE_QUOTE = {".java", ".kt", ".scala", ".swift", ".cs"}
C_FAMILY = C_LIKE | JS_LIKE | TRIPLE_QUOTE | {".go", ".rs", ".sql"}

REGEX_AFTER_CHARS = set("(,=:[!&|?{};+-*%<>~^")
REGEX_AFTER_WORDS = {"return", "typeof", "case", "do", "else", "in", "of", "new",
                     "delete", "void", "throw", "instanceof", "yield", "await"}
IDENT = re.compile(r"[A-Za-z0-9_$]")
WHITESPACE = re.compile(r"\s+")
CPP_RAW_PREFIX = re.compile(r"(?:^|[^A-Za-z0-9_])(?:u8|u|U|L)?R\Z")
RUST_RAW = re.compile(r'b?r(#*)"')
RUST_CHAR = re.compile(r"'(?:[^'\\\n]|\\(?:u\{[0-9A-Fa-f]{1,6}\}|x[0-9A-Fa-f]{2}|.))'")


def _collapse(code):
    # Runs that contain a newline stay a newline: C preprocessor lines end at one.
    return WHITESPACE.sub(lambda m: "\n" if "\n" in m.group() else " ", code)


def strip_c_family(text, ext):
    """Remove comments outside literals, then collapse code whitespace."""
    is_c = ext in C_LIKE
    is_js = ext in JS_LIKE
    is_rs = ext == ".rs"
    is_sql = ext == ".sql"
    pieces = []          # alternating code (to collapse) and literals (verbatim)
    code = []
    recent = ""          # recent significant code, for regex/raw-string context
    i, n = 0, len(text)
    special = re.compile(r"[/\"'`" + ("\\-" if is_sql else "") + ("rb" if is_rs else "") + "]")

    def emit_lit(s):
        nonlocal recent
        pieces.append(_collapse("".join(code)))
        code.clear()
        pieces.append(s)
        recent = '"'     # an operand: a following / is division

    def emit_code(s):
        nonlocal recent
        code.append(s)
        sig = s.strip()
        if sig:
            recent = (recent + s)[-32:]

    while i < n:
        m = special.search(text, i)
        if not m:
            emit_code(text[i:])
            break
        if m.start() > i:
            emit_code(text[i:m.start()])
        i = m.start()
        c = text[i]
        nxt = text[i + 1] if i + 1 < n else ""

        if c == "/" and nxt == "/" or is_sql and c == "-" and nxt == "-":
            j = i
            while True:
                j = text.find("\n", j)
                if j < 0:
                    j = n
                    break
                # C and C++ splice backslash-newline before removing comments.
                if is_c and (text[j - 1] == "\\" or text[j - 2:j] == "\\\r"):
                    j += 1
                    continue
                break
            comment = text[i:j]
            if DIRECTIVES.search(comment):
                emit_lit(comment)
            else:
                code.append(" ")
            i = j
        elif c == "/" and nxt == "*":
            j = text.find("*/", i + 2)
            if j < 0:
                raise Fallback("unterminated block comment")
            comment = text[i:j + 2]
            if DIRECTIVES.search(comment):
                emit_lit(comment)
            else:
                code.append(" ")
            i = j + 2
        elif c == "/" and is_js and _regex_allowed(recent):
            j = i + 1
            in_class = False
            while True:
                if j >= n or text[j] == "\n":
                    raise Fallback("unterminated regex")
                ch = text[j]
                if ch == "\\":
                    j += 2
                    continue
                if ch == "[":
                    in_class = True
                elif ch == "]":
                    in_class = False
                elif ch == "/" and not in_class:
                    break
                j += 1
            j += 1
            while j < n and IDENT.match(text[j]):
                j += 1
            emit_lit(text[i:j])
            i = j
        elif c == "/":
            emit_code(c)
            i += 1
        elif c in "rb" and is_rs and not (i and IDENT.match(text[i - 1])) and RUST_RAW.match(text, i):
            rm = RUST_RAW.match(text, i)
            close = '"' + rm.group(1)
            j = text.find(close, rm.end())
            if j < 0:
                raise Fallback("unterminated raw string")
            emit_lit(text[i:j + len(close)])
            i = j + len(close)
        elif c == '"' and is_c and CPP_RAW_PREFIX.search(recent):
            k = text.find("(", i + 1)
            delim = text[i + 1:k] if k >= 0 else None
            if delim is None or len(delim) > 16 or re.search(r'[\s\\()"]', delim):
                raise Fallback("bad raw string delimiter")
            close = ")" + delim + '"'
            j = text.find(close, k + 1)
            if j < 0:
                raise Fallback("unterminated raw string")
            emit_lit(text[i:j + len(close)])
            i = j + len(close)
        elif c == '"' and ext in TRIPLE_QUOTE and text.startswith('"""', i):
            j = text.find('"""', i + 3)
            if j < 0:
                raise Fallback("unterminated text block")
            emit_lit(text[i:j + 3])
            i = j + 3
        elif c == '"' and ext == ".cs" and recent.endswith("@") and code and code[-1].endswith("@"):
            j = i + 1
            while True:
                j = text.find('"', j)
                if j < 0:
                    raise Fallback("unterminated verbatim string")
                if text.startswith('""', j):
                    j += 2
                    continue
                break
            emit_lit(text[i:j + 1])
            i = j + 1
        elif c == "'" and is_c and _digit_separator(text, i):
            emit_code(c)
            i += 1
        elif c == "'" and is_rs:
            cm = RUST_CHAR.match(text, i)
            if cm:
                emit_lit(cm.group())
                i = cm.end()
            else:            # lifetime or label
                emit_code(c)
                i += 1
        elif c in "\"'":
            j = i + 1
            while True:
                if j >= n:
                    raise Fallback("unterminated string")
                ch = text[j]
                if ch == "\\":
                    j += 2
                    continue
                if ch == "\n":
                    raise Fallback("newline in string")
                if ch == c:
                    break
                j += 1
            emit_lit(text[i:j + 1])
            i = j + 1
        elif c == "`":
            j = _template_end(text, i) if is_js else text.find("`", i + 1)
            if j < 0:
                raise Fallback("unterminated backtick string")
            emit_lit(text[i:j + 1])
            i = j + 1
        else:                # a lone -, r or b that is ordinary code here
            emit_code(c)
            i += 1

    pieces.append(_collapse("".join(code)))
    return "".join(pieces).strip()


def _regex_allowed(recent):
    """Preceding-significant-token heuristic for JS regex literals."""
    s = recent.rstrip()
    if not s:
        return True
    last = s[-1]
    if last in "+-" and len(s) > 1 and s[-2] == last:
        raise Fallback("ambiguous / after ++ or --")
    if last in REGEX_AFTER_CHARS:
        return True
    word = re.search(r"[A-Za-z_$][A-Za-z0-9_$]*\Z", s)
    if word:
        return word.group() in REGEX_AFTER_WORDS
    if IDENT.match(last) or last in ")]\"'`":
        return False
    raise Fallback("cannot classify /")


def _digit_separator(text, i):
    """C++14 digit separators such as 1'000'000."""
    j = i
    while j > 0 and (text[j - 1].isalnum() or text[j - 1] in "_."):
        j -= 1
    return j < i and text[j].isdigit() and i + 1 < len(text) and text[i + 1].isalnum()


def _template_end(text, i):
    """Index of the backtick closing the JS template literal starting at i."""
    j, n, depth = i + 1, len(text), 0
    while j < n:
        ch = text[j]
        if ch == "\\":
            j += 2
            continue
        if depth:
            if ch in "`\"'":
                raise Fallback("quote inside template substitution")
            depth += ch == "{"
            depth -= ch == "}"
        elif ch == "`":
            return j
        elif text.startswith("${", j):
            depth = 1
            j += 1
        j += 1
    return -1


def strip_python(data):
    """Python tokens minus comments and blank-line tokens."""
    out = []
    try:
        for tok in tokenize.tokenize(io.BytesIO(data).readline):
            if tok.type == tokenize.NL:
                continue
            if tok.type == tokenize.COMMENT and not DIRECTIVES.search(tok.string):
                continue
            out.append("%d\x00%s\x01" % (tok.type, tok.string))
    except (tokenize.TokenError, SyntaxError, ValueError) as error:
        raise Fallback(str(error))
    return "".join(out)


YAML_SCALAR_START = set(":-[{,?")
YAML_BLOCK = re.compile(r"(?:^|[\s:\-?])[|>][-+0-9]*\Z")


def strip_hash_comments(text, yaml):
    """Remove # comments from YAML or TOML, string- and block-scalar-aware."""
    lines = text.split("\n")
    out = []              # (line, protected)
    quote = None          # open multi-line quote: '"', "'", '"""' or "'''"
    block_indent = None   # indentation of the line that opened a YAML block scalar
    for line in lines:
        if block_indent is not None:
            stripped = line.lstrip(" ")
            if not stripped.strip() or len(line) - len(stripped) > block_indent:
                out.append((line, True))
                continue
            block_indent = None
        protected = quote is not None
        kept, directive, i, n = [], "", 0, len(line)
        while i < n:
            c = line[i]
            if quote:
                if quote in ('"""', "'''") and line.startswith(quote, i):
                    kept.append(quote)
                    i += 3
                    quote = None
                    continue
                if c == "\\" and quote[0] == '"':
                    kept.append(line[i:i + 2])
                    i += 2
                    continue
                if yaml and quote == "'" and line.startswith("''", i):
                    kept.append("''")
                    i += 2
                    continue
                kept.append(c)
                i += 1
                if len(quote) == 1 and c == quote:
                    quote = None
                continue
            if c == "#" and (not yaml or i == 0 or line[i - 1] in " \t"):
                if DIRECTIVES.search(line[i:]):
                    directive = line[i:]
                break
            if c in "\"'":
                before = "".join(kept).rstrip()
                if not yaml:
                    quote = c * 3 if line.startswith(c * 3, i) else c
                    kept.append(quote)
                    i += len(quote)
                    continue
                if not before.strip() or before[-1] in YAML_SCALAR_START:
                    quote = c
            kept.append(c)
            i += 1
        if quote is not None:
            if not yaml and len(quote) == 1:
                raise Fallback("unterminated TOML string")
            protected = True
        text_kept = "".join(kept)
        out.append((text_kept + directive, protected))
        if yaml and quote is None and YAML_BLOCK.search(text_kept.rstrip()):
            block_indent = len(line) - len(line.lstrip(" "))
    if quote is not None:
        raise Fallback("unterminated string")
    result = []
    for line, protected in out:
        if protected:
            result.append(line)
        elif line.strip():
            result.append(line.rstrip())
    return "\n".join(result)


def normalize(path, data, strip_comments):
    """Bytes whose equality means 'no semantic change' for this file."""
    if not strip_comments:
        return data
    ext = os.path.splitext(path)[1]
    try:
        if ext == ".py":
            return b"py\x00" + strip_python(data).encode("utf-8")
        text = data.decode("utf-8")
        if ext in C_FAMILY:
            return b"c\x00" + strip_c_family(text, ext).encode("utf-8")
        if ext in (".yml", ".yaml"):
            return b"yaml\x00" + strip_hash_comments(text, yaml=True).encode("utf-8")
        if ext == ".toml":
            return b"toml\x00" + strip_hash_comments(text, yaml=False).encode("utf-8")
    except Exception:  # noqa: BLE001 - any doubt means compare raw bytes
        return data
    return data


# --------------------------------------------------------------------------
# Git

def git(*args, data=None):
    result = subprocess.run(("git",) + args, input=data, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, check=True)
    return result.stdout


class Blobs:
    """Reads blobs through one `git cat-file --batch` process."""

    def __init__(self):
        self.proc = subprocess.Popen(["git", "cat-file", "--batch"], stdin=subprocess.PIPE,
                                     stdout=subprocess.PIPE)

    def read(self, sha):
        self.proc.stdin.write(sha.encode() + b"\n")
        self.proc.stdin.flush()
        header = self.proc.stdout.readline().split()
        if len(header) < 3 or header[1] == b"missing":
            raise RuntimeError("missing object %s" % sha)
        size = int(header[2])
        body = self.proc.stdout.read(size)
        self.proc.stdout.read(1)
        return body

    def close(self):
        self.proc.stdin.close()
        self.proc.stdout.close()
        self.proc.wait()


def changed_entries(base, head):
    """(path, old mode, new mode, old sha, new sha) for base...head."""
    raw = git("diff", "--raw", "-z", "--no-renames", "--no-abbrev", "%s...%s" % (base, head))
    fields = raw.split(b"\0")
    entries = []
    for k in range(0, len(fields) - 1, 2):
        meta = fields[k].decode().lstrip(":").split()
        entries.append((fields[k + 1].decode("utf-8", "surrogateescape"),
                        meta[0], meta[1], meta[2], meta[3]))
    return entries


def tree_entries(rev):
    """(path, mode, type, sha) for every entry in rev's tree."""
    raw = git("ls-tree", "-r", "-z", "--full-tree", rev)
    entries = []
    for record in raw.split(b"\0"):
        if not record:
            continue
        meta, path = record.split(b"\t", 1)
        mode, kind, sha = meta.decode().split()
        entries.append((path.decode("utf-8", "surrogateescape"), mode, kind, sha))
    return entries


# --------------------------------------------------------------------------
# Areas

class Config:
    def __init__(self, raw):
        self.raw = raw
        spec = json.loads(raw)
        self.global_ = Matcher(spec.get("global", []))
        self.ignore = Matcher(spec.get("ignore", []))
        self.areas = {}
        for name, area in spec["areas"].items():
            if not re.match(r"[a-z][a-z0-9_]*\Z", name):
                raise ValueError("area names must be lowercase identifiers: %r" % name)
            self.areas[name] = (Matcher(area["paths"]), bool(area.get("strip_comments", True)))

    @classmethod
    def load(cls, path=CONFIG):
        with open(path, "rb") as handle:
            return cls(handle.read())

    def in_footprint(self, area, path):
        """Global files are in every footprint; ignored files in none."""
        if self.global_(path):
            return True
        return not self.ignore(path) and self.areas[area][0](path)


def detect_changes(config, base, head, blobs):
    """Map area -> (affected, reason)."""
    everything = lambda reason: {a: (True, reason) for a in config.areas}
    if not base or set(base) == {"0"}:
        return everything("no base commit")
    try:
        entries = changed_entries(base, head)
    except (subprocess.CalledProcessError, OSError, ValueError, IndexError):
        return everything("git diff failed")
    for path, *_ in entries:
        if config.global_(path):
            return everything("global: " + path)
    for path, old_mode, new_mode, _, _ in entries:
        if config.ignore(path):
            # Docs may be installed or asserted to exist: only edits are free.
            if "000000" in (old_mode, new_mode):
                return everything("added/removed: " + path)
        elif not any(m(path) for m, _ in config.areas.values()):
            return everything("unmapped: " + path)

    result = {}
    for area, (_, strip) in config.areas.items():
        result[area] = (False, "no relevant changes")
        for path, old_mode, new_mode, old_sha, new_sha in entries:
            if not config.in_footprint(area, path):
                continue
            try:
                if old_mode != new_mode or _content(blobs, path, old_sha, strip) != \
                        _content(blobs, path, new_sha, strip):
                    result[area] = (True, "changed: " + path)
                    break
            except Exception:  # noqa: BLE001
                result[area] = (True, "unreadable: " + path)
                break
        else:
            if any(config.in_footprint(area, e[0]) for e in entries):
                result[area] = (False, "comment/whitespace-only changes")
    return result


def _content(blobs, path, sha, strip):
    if set(sha) == {"0"}:
        return None
    return normalize(path, blobs.read(sha), strip)


def area_keys(config, head, blobs, extra_files, toolchain):
    """Map area -> SHA-256 fingerprint of its normalized inputs at head."""
    entries = sorted(tree_entries(head))
    extras = b"".join(b"%d\0" % k + _read(p) for k, p in enumerate(extra_files))
    digest_cache = {}
    keys = {}
    for area, (_, strip) in config.areas.items():
        h = hashlib.sha256()
        for part in (SCRIPT_VERSION, area, sys.version, toolchain):
            h.update(part.encode() + b"\0")
        h.update(extras)
        for path, mode, kind, sha in entries:
            if not config.in_footprint(area, path):
                if config.ignore(path):
                    h.update(b"ignored\0" + path.encode("utf-8", "surrogateescape") + b"\n")
                continue
            # Global files are compared raw, as in change detection.
            use_strip = strip and not config.global_(path) and kind == "blob" and mode != "120000"
            cache_key = (sha, use_strip, os.path.splitext(path)[1])
            if cache_key not in digest_cache:
                if kind != "blob":
                    digest_cache[cache_key] = sha.encode()
                else:
                    content = blobs.read(sha)
                    digest_cache[cache_key] = hashlib.sha256(
                        normalize(path, content, use_strip)).hexdigest().encode()
            h.update(path.encode("utf-8", "surrogateescape") + b"\0" + mode.encode() + b"\0")
            h.update(digest_cache[cache_key] + b"\n")
        keys[area] = h.hexdigest()
    return keys


def _read(path):
    try:
        with open(path, "rb") as handle:
            return handle.read() + b"\0"
    except OSError:
        return b"<missing>\0"


# --------------------------------------------------------------------------
# Commands

def _clean(text):
    return re.sub(r"[^\w./:+\- ]", "_", text)[:120]


def _summary(lines):
    path = os.environ.get("GITHUB_STEP_SUMMARY")
    if path:
        with open(path, "a") as handle:
            handle.write("\n".join(lines) + "\n\n")


def fingerprint(args):
    config = Config.load(args.config)
    toolchain = os.environ.get("TOOLCHAIN", "")
    blobs = Blobs()
    try:
        try:
            changes = detect_changes(config, args.base, args.head, blobs)
        except Exception as error:  # noqa: BLE001
            changes = {a: (True, "detection error: %s" % error) for a in config.areas}
        try:
            keys = area_keys(config, args.head, blobs,
                             [args.config, os.path.abspath(__file__), args.workflow], toolchain)
        except Exception as error:  # noqa: BLE001
            # A random key never matches a recorded pass marker.
            print("::warning::fingerprint failed (%s); no pass markers will match" % error,
                  file=sys.stderr)
            keys = {a: "unmatched-" + uuid.uuid4().hex for a in config.areas}
    finally:
        blobs.close()

    rows = ["### Change detection", "",
            "| Area | Affected | Reason | Key |", "|---|---|---|---|"]
    for area in config.areas:
        affected, reason = changes[area]
        print("%s_affected=%s" % (area, "true" if affected else "false"))
        print("%s_key=%s" % (area, keys[area]))
        print("%s_reason=%s" % (area, _clean(reason)))
        rows.append("| %s | %s | %s | `%s` |" % (area, "yes" if affected else "no",
                                                 _clean(reason), keys[area][:12]))
    _summary(rows)


def decide(args):
    """area_run = affected and no pass marker for this exact key."""
    config = Config.load(args.config)
    force = os.environ.get("FORCE", "").lower() == "true"
    rows = ["### Decision", "", "| Area | Affected | Pass marker | Runs |", "|---|---|---|---|"]
    for area in config.areas:
        env = area.upper()
        affected = os.environ.get("AFFECTED_" + env, "") != "false"
        key = os.environ.get("KEY_" + env, "")
        hit = os.environ.get("HIT_" + env, "") == "true" and re.match(r"[0-9a-f]{64}\Z", key)
        run = force or (affected and not hit)
        print("%s_run=%s" % (area, "true" if run else "false"))
        rows.append("| %s | %s | %s | %s |" % (
            area, "yes" if affected else "no", "hit" if hit else "-",
            ("**yes** (forced)" if force else "**yes**") if run else "skipped"))
    _summary(rows)


def check_workflow(config, text):
    """Problems that would let ci.yml drift from areas.json."""
    problems = []
    for area in config.areas:
        env = area.upper()
        expected = [
            ("changes output %s_run" % area,
             r"%s_run: \$\{\{ steps\.decide\.outputs\.%s_run \}\}" % (area, area)),
            ("changes output %s_key" % area,
             r"%s_key: \$\{\{ steps\.fp\.outputs\.%s_key \}\}" % (area, area)),
            ("marker lookup step pass_%s" % area,
             r"id: pass_%s\b[\s\S]*?key: pass-%s-\$\{\{ steps\.fp\.outputs\.%s_key \}\}"
             % (area, area, area)),
            ("decide env for %s" % area,
             r"AFFECTED_%s: \$\{\{ steps\.fp\.outputs\.%s_affected \}\}[\s\S]*?"
             r"KEY_%s: \$\{\{ steps\.fp\.outputs\.%s_key \}\}[\s\S]*?"
             r"HIT_%s: \$\{\{ steps\.pass_%s\.outputs\.cache-hit \}\}"
             % (env, area, env, area, env, area)),
            ("job gate for %s" % area,
             r"needs\.changes\.outputs\.%s_run == 'true'" % area),
            ("marker save for %s" % area,
             r"key: pass-%s-\$\{\{ needs\.changes\.outputs\.%s_key \}\}" % (area, area)),
            ("status dependency on %s" % area,
             r"(?m)^  status:[\s\S]*?^    - %s$" % area),
        ]
        for label, pattern in expected:
            if not re.search(pattern, text):
                problems.append("missing " + label)
    gated = set(re.findall(r"needs\.changes\.outputs\.([a-z0-9_]+)_run == 'true'", text))
    for extra in sorted(gated - set(config.areas)):
        problems.append("job gated on unknown area %s" % extra)
    return problems


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--base", default="")
    parser.add_argument("--head", default="HEAD")
    parser.add_argument("--config", default=CONFIG)
    parser.add_argument("--workflow", default=WORKFLOW)
    parser.add_argument("--decide", action="store_true")
    parser.add_argument("--check-workflow", metavar="PATH")
    args = parser.parse_args(argv)
    if args.check_workflow:
        with open(args.check_workflow) as handle:
            problems = check_workflow(Config.load(args.config), handle.read())
        for problem in problems:
            print("::error file=%s::%s" % (args.check_workflow, problem))
        return 1 if problems else 0
    if args.decide:
        decide(args)
    else:
        fingerprint(args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
