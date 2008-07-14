/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2005-2008 by Lukas Gamper <mistral@student.ethz.ch>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* This software is part of the ALPS libraries, published under the ALPS
* Library License; you can use, redistribute it and/or modify it under
* the terms of the license, either version 1 or (at your option) any later
* version.
* 
* You should have received a copy of the ALPS Library License along with
* the ALPS Libraries; see the file LICENSE.txt. If not, the license is also
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

/* $Id$ */

#ifndef _NODE_H_
#define _NODE_H_

#include <string>
#include <list>
#include <map>
#include <algorithm>
#include <iostream>

class Node {
        std::string mName;
        Node *mParent;
        std::map<std::string, std::string> mAttr;
        std::list<std::string> mText;
        std::list<Node> mElem;

        public:
                Node(): mName(""), mParent(NULL) {}
                Node(std::string inName, Node *inParent): mName(inName), mParent(inParent) {}

                void setParent(Node *inParent) { mParent = inParent; }
                Node *getParent() { return mParent; }

                void setName(std::string inName) { mName = inName; }
                std::string getName() { return mName; }

                void addElement(Node inChild);
                void addElement(std::string inName);
                std::list<Node> getElements() { return mElem; }
                Node *getLastElement() { return &mElem.back(); }

                void addText(std::string inName) { mText.push_back(inName); }
                std::list<std::string> getText() const { return mText; }

                void addAttribute(std::string inName, std::string inValue);
                std::string getAttribute(std::string inName) const;

                std::string string();
                Node getElement(std::string nodeName);
                std::list<Node> nodeTest(std::string nodeName);
};
#endif //_NODE_H_
