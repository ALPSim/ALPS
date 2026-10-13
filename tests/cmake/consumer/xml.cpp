// SPDX-License-Identifier: MIT
#include <alps/parser/parser.h>
#include <alps/xml.h>
#include <sstream>
#include <stdexcept>

int main() {
    // Exercise writer, attribute parsing and the handler-based reader through
    // the installed XML component, without loading simulation or archive code.
    std::ostringstream output;
    {
        alps::oxstream xml(output);
        xml << alps::start_tag("VALUE") << alps::attribute("name", "energy")
            << 42 << alps::end_tag("VALUE");
    }
    std::istringstream tags(output.str());
    const auto tag = alps::parse_tag(tags);
    if (tag.name != "VALUE" || tag.attributes["name"] != "energy")
        throw std::runtime_error("XML component attribute contract failed");
    int value = 0;
    alps::SimpleXMLHandler<int> handler("VALUE", value);
    alps::XMLParser parser(handler);
    std::istringstream input(output.str());
    parser.parse(input);
    if (value != 42)
        throw std::runtime_error("XML component round-trip contract failed");
}
