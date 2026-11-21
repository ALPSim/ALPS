// Copyright (C) 2026 ALPS Collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT

#ifndef ALPS_DMFT_TEXT_TABLE_H
#define ALPS_DMFT_TEXT_TABLE_H

/// @file text_table.h
/// @brief reads a text table of numbers and checks that it is complete
///
/// The Green's function and hybridization function text readers read a fixed
/// number of rows and stopped there: a file with more rows was cut, one with
/// fewer left the rest at 0 or uninitialized, and a value that was not a
/// number ended the reading without a message.

#include <complex>  // read_table<std::complex<double> > reads "(re,im)"
#include <istream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace text_table {

/// Reads a table of `rows` rows of one label column (not used) and `columns`
/// values, and returns the values row by row.
template <class T>
std::vector<T> read_table(std::istream& in, std::size_t rows, std::size_t columns,
                          std::string const& what, std::string const& file)
{
  // the label column is real also in a complex table
  std::vector<T> values;
  values.reserve(rows * columns);
  std::size_t count = 0;
  std::string word;
  std::vector<std::string> words;
  while (in >> word)
    words.push_back(word);
  if (words.size() != rows * (columns + 1)) {
    std::ostringstream msg;
    msg << what << " " << file << ": has " << words.size() << " entries; it needs " << rows
        << " rows of " << columns + 1 << " (a label and " << columns << " values), " << rows * (columns + 1) << ".";
    throw std::runtime_error(msg.str());
  }
  for (std::size_t r = 0; r < rows; ++r)
    for (std::size_t c = 0; c <= columns; ++c) {
      std::istringstream entry(words[count++]);
      if (c == 0) {
        double label;
        if (!(entry >> label) || !(entry >> std::ws).eof())
          throw std::runtime_error(what + " " + file + ": entry " + std::to_string(count) + " (\""
                                   + words[count - 1] + "\") is not a number.");
        continue;
      }
      T x;
      if (!(entry >> x) || !(entry >> std::ws).eof())
        throw std::runtime_error(what + " " + file + ": entry " + std::to_string(count) + " (\""
                                 + words[count - 1] + "\") is not a number.");
      values.push_back(x);
    }
  return values;
}

}  // namespace text_table

#endif
