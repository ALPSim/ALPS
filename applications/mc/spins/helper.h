/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C)  by Andreas Streich 
*
* This software is part of the ALPS Applications, published under the ALPS
* Application License; you can use, redistribute it and/or modify it under
* the terms of the license, either version 1 or (at your option) any later
* version.
* 
* You should have received a copy of the ALPS Application License along with
* the ALPS Applications; see the file LICENSE.txt. If not, the license is also
* available from http://alps.comp-phys.org/.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
* FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
* SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
* FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
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
