// SPDX-License-Identifier: MIT
#include <alps/cli/options.hpp>

#include <algorithm>
#include <ostream>

namespace alps {
namespace cli {

namespace {
// The error texts follow Boost.program_options, which these programs used before.
std::string dashed(const std::string& name) { return "'--" + name + "'"; }
}

options::options(const std::string& caption) : caption_(caption) {}

options& options::section(const std::string& caption) {
  caption_ = caption;
  return *this;
}

options& options::add(const std::string& names, const std::string& help,
                      std::function<bool(const std::string&)> assign) {
  const std::string::size_type comma = names.find(',');
  option opt;
  opt.name = names.substr(0, comma);
  opt.letter = comma == std::string::npos ? '\0' : names.at(comma + 1);
  opt.help = help;
  opt.caption = caption_;
  opt.takes_value = static_cast<bool>(assign);
  opt.repeatable = opt.hidden = false;
  opt.assign = std::move(assign);
  opt.seen = 0;
  options_.push_back(std::move(opt));
  return *this;
}

options& options::flag(const std::string& names, const std::string& help) {
  return add(names, help, nullptr);
}

options& options::list(const std::string& names, const std::string& help) {
  add(names, help, [](const std::string&) { return true; });
  options_.back().repeatable = true;
  return *this;
}

options& options::hidden() {
  options_.back().hidden = true;
  return *this;
}

options& options::fallback_text(const std::string& text) {
  options_.back().fallback = text;
  return *this;
}

options& options::positional(const std::string& name, int max_count) {
  find(name);
  positional_.emplace_back(name, max_count);
  return *this;
}

const options::option& options::find(const std::string& name) const {
  for (const option& opt : options_)
    if (opt.name == name) return opt;
  throw error("option " + dashed(name) + " is not declared");
}

options::option& options::find(const std::string& name) {
  return const_cast<option&>(static_cast<const options&>(*this).find(name));
}

std::size_t options::count(const std::string& name) const {
  // An option that a build does not declare (--mpi, --threaded) was not given.
  for (const option& opt : options_)
    if (opt.name == name) return opt.seen;
  return 0;
}

const std::vector<std::string>& options::values(const std::string& name) const { return find(name).given; }

options::option& options::match_long(const std::string& name) {
  option* match = nullptr;
  std::vector<std::string> names;
  for (option& opt : options_) {
    if (opt.name == name) return opt;
    if (!name.empty() && opt.name.compare(0, name.size(), name) == 0) {
      match = &opt;
      names.push_back(dashed(opt.name));
    }
  }
  if (names.empty()) throw error("unrecognised option " + dashed(name));
  if (names.size() == 1) return *match;
  std::sort(names.begin(), names.end());
  std::string message = "option " + dashed(name) + " is ambiguous and matches ";
  for (std::size_t i = 0; i < names.size(); ++i)
    message += (i == 0 ? "" : i + 1 == names.size() ? ", and " : ", ") + names[i];
  throw error(message);
}

options::option* options::match_short(char letter) {
  for (option& opt : options_)
    if (opt.letter == letter) return &opt;
  return nullptr;
}

void options::store(option& opt, const std::string* text) {
  if (opt.seen && !opt.repeatable)
    throw error("option " + dashed(opt.name) + " cannot be specified more than once");
  ++opt.seen;
  if (text) opt.given.push_back(*text);
}

void options::parse(int argc, const char* const* argv) {
  const std::vector<std::string> args(argv + (argc > 0 ? 1 : 0), argv + argc);
  bool options_ended = false;
  std::size_t next_positional = 0;
  for (std::size_t i = 0; i < args.size(); ++i) {
    const std::string& arg = args[i];
    // The value of an option is the next argument, whatever it looks like,
    // except a short option on its own ("-T -c" is a missing value).
    auto following = [&](const option& opt) -> const std::string* {
      const bool is_short = i + 1 < args.size() && args[i + 1].size() == 2 && args[i + 1][0] == '-'
                            && match_short(args[i + 1][1]);
      if (i + 1 == args.size() || is_short)
        throw error("the required argument for option " + dashed(opt.name) + " is missing");
      return &args[++i];
    };
    if (options_ended || arg.size() < 2 || arg[0] != '-') {
      std::size_t slot = 0, before = 0;
      while (slot < positional_.size() && positional_[slot].second >= 0
             && next_positional >= before + positional_[slot].second)
        before += positional_[slot++].second;
      if (slot == positional_.size())
        throw error("too many positional options have been specified on the command line");
      ++next_positional;
      store(find(positional_[slot].first), &arg);
    } else if (arg == "--") {
      options_ended = true;
    } else if (arg[1] == '-') {
      const std::string::size_type equals = arg.find('=');
      option& opt = match_long(arg.substr(2, equals == std::string::npos ? equals : equals - 2));
      if (!opt.takes_value) {
        if (equals != std::string::npos)
          throw error("option " + dashed(opt.name) + " does not take any arguments");
        store(opt, nullptr);
      } else if (equals == std::string::npos) {
        store(opt, following(opt));
      } else {
        if (equals + 1 == arg.size())
          throw error("the argument for option " + dashed(opt.name) + " should follow immediately after the equal sign");
        const std::string text = arg.substr(equals + 1);
        store(opt, &text);
      }
    } else {
      for (std::size_t j = 1; j < arg.size(); ++j) {
        option* opt = match_short(arg[j]);
        if (!opt) throw error("unrecognised option '" + arg + "'");
        if (!opt->takes_value) {
          store(*opt, nullptr);
          continue;
        }
        const std::string rest = arg.substr(j + 1);
        store(*opt, rest.empty() ? following(*opt) : &rest);
        break;
      }
    }
  }
  // Values are converted only once the whole command line is well formed.
  for (option& opt : options_)
    for (const std::string& text : opt.given)
      if (!opt.assign(text))
        throw error("the argument ('" + text + "') for option " + dashed(opt.name) + " is invalid");
}

void options::print(std::ostream& os) const {
  std::vector<std::string> entries;
  std::size_t width = 21;  // the narrowest column Boost.program_options used
  for (const option& opt : options_) {
    std::string entry = opt.letter ? std::string("-") + opt.letter + " [ --" + opt.name + " ]" : "--" + opt.name;
    if (opt.takes_value) entry += " arg";
    if (!opt.fallback.empty()) entry += " (=" + opt.fallback + ")";
    if (!opt.hidden) width = (std::max)(width, entry.size() + (opt.takes_value ? 0 : 1));
    entries.push_back(entry);
  }
  const std::size_t column = width + 3, line_length = 79;
  const std::string* caption = nullptr;
  for (std::size_t i = 0; i < options_.size(); ++i) {
    const option& opt = options_[i];
    if (opt.hidden) continue;
    if (!caption || *caption != opt.caption) {
      if (caption) os << "\n";
      if (!opt.caption.empty()) os << opt.caption << ":\n";
      caption = &opt.caption;
    }
    os << "  " << entries[i] << std::string(column - 2 - entries[i].size(), ' ');
    // Break the description between words where it would pass the line length.
    std::istringstream words(opt.help);
    std::string word;
    std::size_t used = column;
    for (bool first = true; words >> word; first = false) {
      if (!first && used + 1 + word.size() > line_length) {
        os << "\n" << std::string(column, ' ');
        used = column;
      } else if (!first) {
        os << ' ';
        ++used;
      }
      os << word;
      used += word.size();
    }
    os << "\n";
  }
}

} // namespace cli
} // namespace alps
