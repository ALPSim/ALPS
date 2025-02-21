/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C)  by Andreas Streich 
*
* Permission is hereby granted, free of charge, to any person obtaining
* a copy of this software and associated documentation files (the “Software”),
* to deal in the Software without restriction, including without limitation
* the rights to use, copy, modify, merge, publish, distribute, sublicense,
* and/or sell copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included
* in all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS
* OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
* DEALINGS IN THE SOFTWARE.
*
*****************************************************************************/

#ifndef ALPS_APPLICATIONS_MC_SPIN_HELPER_H_
#define ALPS_APPLICATIONS_MC_SPIN_HELPER_H_

#include <algorithm>
#include <iostream>
#include <alps/parameter.h>

/* astreich, 06/14 */
template<class T>
class Helper {
public:
  Helper() {};
  Helper(std::string& str) { 
    myString = str;
    myPos = 0;
    paramsSet = false;
  }
  Helper(std::string& str,alps::Parameters& params) {
    myString = str;
    myPos = 0;
    myParams = params;
    paramsSet = true;
  }
 
  void setString(const std::string& str) {
    myString = str;
    myPos = 0;
    paramsSet = false;
  }
 
  int elemcount(std::string& str) {
    int counter = 0;
    int pos = 0;
    if (str[0] == '\0') return 0;
    while (true) {
      while (str[pos] == ' ') pos++;
      if (str[pos] == '\0') return counter;
      counter++;
      while ((str[pos] !='\0') && (str[pos] != ' ')) pos++;
    }
  }

  int elemcount() {
    int counter = 0;
    int pos = 0;
    if (myString[0] == '\0') return 0;
    while (true) {
      while (myString[pos] == ' ') pos++;
      if (myString[pos] == '\0') return counter;
      counter++;
      while ((myString[pos] !='\0') && (myString[pos] != ' ')) pos++;
    }
  }

  T getNextElement() {
    std::string tmp;
    tmp.clear();
    while (myString[myPos] == ' ') { myPos++; }
    while ((myString[myPos] != ' ') && (myString[myPos] != '\0')) {
      tmp.append(1,myString[myPos]);
      myPos++;
    }
    if (tmp.length()==0) {
      std::cerr << "no more arguments left\n";
      boost::throw_exception(std::runtime_error("too few input values"));
    }
    if (paramsSet) return alps::evaluate<T>(tmp, myParams);
    else return alps::evaluate<T>(tmp);
  }

  std::string getNextString() {
    std::string res;
    res.clear();
    while (myString[myPos] == ' ') myPos++;
    while ((myString[myPos] != ' ') && (myString[myPos] != '\0')) {
      res.append(1,myString[myPos]);
      myPos++;
    }
    return res;
  }
private:
  std::string myString;
  int myPos;
  alps::Parameters myParams;
  bool paramsSet;
};

#endif
