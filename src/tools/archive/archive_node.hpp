/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2005-2008 by Lukas Gamper <mistral@student.ethz.ch>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
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
