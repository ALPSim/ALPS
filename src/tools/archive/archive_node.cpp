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

#include <stdexcept>

#include "archive_node.hpp"

void Node::addElement(Node inChild) {
        inChild.mParent = this;
        mElem.push_back(inChild);
}

void Node::addElement(std::string inName) {
        std::transform (inName.begin(), inName.end(), inName.begin(), tolower);
        addElement(Node(inName, this));
}

void Node::addAttribute(std::string inName, std::string inValue) {
        std::transform (inName.begin(), inName.end(), inName.begin(), tolower);
        mAttr[inName] = inValue;
}

std::string Node::getAttribute(std::string inName) const {
  std::transform (inName.begin(), inName.end(), inName.begin(), tolower);
  std::map<std::string, std::string>::const_iterator itr = mAttr.find(inName);
  if (itr == mAttr.end())
    return "";
  else
    return itr->second;
}

std::string Node::string() {
        std::string value;
        for (std::list<std::string>::iterator it = mText.begin(); it != mText.end(); it++)
                value += *it;
        return value;
}

std::list<Node> Node::nodeTest(std::string nodeName) {
        std::transform (nodeName.begin(), nodeName.end(), nodeName.begin(), tolower);
        std::list<Node> context;
        for (std::list<Node>::iterator it = mElem.begin(); it != mElem.end(); it++)
                if (it->getName() == nodeName)
                        context.push_back(*it);
        return context;
}

Node Node::getElement(std::string nodeName) {
        std::transform (nodeName.begin(), nodeName.end(), nodeName.begin(), tolower);
        for (std::list<Node>::iterator it = mElem.begin(); it != mElem.end(); it++)
                if (it->getName() == nodeName)
                        return *it;
        return Node();
}
