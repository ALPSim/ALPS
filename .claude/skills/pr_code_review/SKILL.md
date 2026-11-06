---
name: "pr_code_review"
description: "Organize a PR's changed files into a categorized HTML review page. Run with /pr_code_review <PR number or URL>."
disable-model-invocation: true
---

# PR review: group the changes for easier review

Input: a GitHub PR, given as a number or URL.

Goal: show every changed file in the PR exactly once, grouped into categories so it is easier to review, and open the result in a new browser tab. You are organizing the changes, not judging them. The reviewer decides what matters.

## 1. Get the PR

```bash
gh pr view <PR> --json number,title,body,url,headRefOid
gh pr diff <PR> > pr.diff
```

Save `pr.diff` and any scripts outside the repo. Read the title and description first.

**Stacked PRs.** If the PR is stacked on other open PRs (the description says "stacked on #N", or the commits of another open PR are in its branch), GitHub's diff against the base branch also contains every PR below it in the stack. Review only this PR's own changes: diff its head against the head of the PR directly below it, which already includes the ones further down.

```bash
gh pr view <lower PR> --json headRefName,headRefOid      # the PR directly below
git fetch origin <lower head branch> <this head branch>
git diff <lower head sha> <this head sha> > pr.diff
git diff --numstat <lower head sha> <this head sha>      # expected files and line counts
```

Use this `pr.diff` everywhere below. Lines removed in it that the lower PR added mean this branch is missing a later change from the lower PR, not that this PR removes them: check with `git merge-base --is-ancestor`, and list such files in chat instead of reviewing them as this PR's changes. On the page, say which PR it is stacked on and that its changes are left out.

## 2. Read every change

Go through the whole diff, hunk by hunk. Don't skim, and don't sample.

## 3. Group the files into categories

The unit is the **file**: every changed file goes into exactly one category, with all of its hunks kept together. Don't split a file's hunks across categories.

Work out the categories from this PR. There is no fixed list. Group files by how similar their changes are and how they depend on each other, and name each category after what its changes actually do. For example:
- a mechanical change repeated across many files (such as replacing `parms["KEY"]` map lookups with fields of a struct) is one category, with the file that defines the struct or helper first, followed by every file that switches to it;
- a change to one program's behavior (such as how a tool reads its input) is its own category.

Rules:
- **One-to-one.** Every file in the diff appears in exactly one category. Nothing is left out and nothing appears twice.
- **Stackable categories.** Shape the categories so each one could become its own smaller PR, stacked on the ones it depends on. If a category would depend on another only through one or two of its files, move those files or split the category, so that the dependencies between categories stay few and clear.
- **Find the dependencies from the code.** Category B depends on category A when B's files use something A adds or changes: a new header, function, type or signature, or a behavior that B's tests check. Check this in the diff (includes, symbol uses, test expectations), and check the PR branch when a dependency runs through code outside the diff, such as an include that is unchanged. Don't guess from file names.
- **Number the categories** in dependency order, so each one comes after everything it depends on.
- **Two categories go last**, after every other category including miscellaneous, so the less trivial changes are reviewed first:
  1. **Simple swaps**: files whose only change is swapping one type, name or spelling for another, such as `int` to `std::size_t`, with no other change in behavior. If another category depends on a swap (for example, a function signature that other files now call with the new type), it is not a simple swap: keep it in dependency order with the change that needs it.
  2. **Tests and benchmarks** (the very last category): test files, test inputs and reference outputs, test CMake files, and benchmark scripts and results. A test file that also changes what is built or installed still goes here, and its note names that part.
  Leave either category out if the PR has no such files.
- **Order within a category.** Put the files that define something (headers, new helpers) before the files that use it.
- **A file with mixed changes** goes in the category of its main change. Its note names the other parts and where they are in the file.
- **Leftovers are fine.** Files that don't fit any category cleanly go together in a miscellaneous category, placed just before the simple swaps and tests and benchmarks.
- Show every file in full, even in a repetitive category. Don't call any change unimportant, safe, or skippable.

### Descriptions and notes

The reader has limited attention and the diff is right below your text. Write
what the diff does not make obvious, in plain sentences a person would say out
loud, not an inventory of everything that changed.

- **Category description:** one sentence: what the change does and why.
- **File note:** one or two short sentences. Start with the main change in this
  file. Then name only what a reader could miss in the diff: a part that differs
  from the rest of its category, something removed without replacement, or
  where code moved to or from.
- Don't list every section, function or line touched; the diff shows that.
- **Depends on / Needed by / Shared files:** one short line per entry: the
  category, the thing it uses, and the files that use it.

If a note needs more than two sentences, the file likely belongs in its own
category, or the extra detail belongs in the diff, not the note.

## 4. Check coverage

Count the files in `pr.diff` and the files across all categories. The two numbers must match, with no duplicates. Also check that the page shows the same number of added and removed lines as the PR (for a stacked PR, as `git diff --numstat` against the lower PR's head). If anything doesn't match, fix the grouping.

## 5. Build the page

Make one self-contained HTML file containing:
- the PR title, number, and a link to the PR
- a numbered list of the categories with their file counts and what each depends on, each linking to its section
- for each category, before its files:
  - **Depends on:** the categories it needs, and for each one what it uses from it (the header, function or behavior) and which of its files use it
  - **Needed by:** the categories that depend on it, with the files in them that use it
  - **Shared files:** any file in it that holds changes for several categories (such as a CMakeLists.txt, a changelog, or a test file covering several parts), and what would have to be split out of it if the categories became separate PRs
  - a category that depends on nothing and that nothing depends on says so
- for each category: its name, its description, and each file's full diff. Each file gets its path, its note, a link to that file in the PR's "Files changed" tab, and its hunks shown in order, each with old and new line numbers, red for removed lines and green for added lines
- new files collapsed by default, expandable with a click
- a "reviewed" checkbox on each file, saved in the browser and keyed to the PR's head commit so checks reset when the PR changes

Style it like GitHub's diff view (plain, monospace code, light background).

Generate the diffs with a small script that reads `pr.diff`, not by copying code by hand, so the code shown is exactly what's in the PR. Only the category assignments and notes come from you. The script should refuse to build the page if the coverage check in step 4 fails.

## 6. Open it

Save it as `pr-<number>-review.html` outside the repo and open it in the browser (`open` on macOS, `xdg-open` on Linux, `start` on Windows).

In chat, list the numbered categories with their file counts and what each depends on, and give the file path. Don't post anything to GitHub.
