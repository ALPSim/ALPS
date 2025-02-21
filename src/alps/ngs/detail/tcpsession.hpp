/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2011 by Lukas Gamper <gamperl@gmail.com>                   *
 *                                                                                 *
 * Permission is hereby granted, free of charge, to any person obtaining           *
 * a copy of this software and associated documentation files (the “Software”),    *
 * to deal in the Software without restriction, including without limitation       *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,        *
 * and/or sell copies of the Software, and to permit persons to whom the           *
 * Software is furnished to do so, subject to the following conditions:            *
 *                                                                                 *
 * The above copyright notice and this permission notice shall be included         *
 * in all copies or substantial portions of the Software.                          *
 *                                                                                 *
 * THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS         *
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,     *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE     *
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER          *
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING         *
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER             *
 * DEALINGS IN THE SOFTWARE.                                                       *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef ALPS_NGS_DETAIL_TCPSESSION_HPP
#define ALPS_NGS_DETAIL_TCPSESSION_HPP

#include <boost/bind.hpp>
#include <boost/asio.hpp>
#include <boost/lexical_cast.hpp>

#include <cstdlib>

namespace alps {

    namespace detail {

        class tcpsession {

            public:

                typedef std::map<std::string, boost::function<std::string()> > actions_type;

                tcpsession(actions_type & a, boost::asio::io_service & io_service)
                    : actions(a)
                    , socket(io_service)
                {}

                boost::asio::ip::tcp::socket & get_socket() {
                    return socket;
                }

                void start(boost::system::error_code const & error = boost::system::error_code()) {
                    if (!error)
                        socket.async_read_some(
                              boost::asio::buffer(data, 1028)
                            , boost::bind(
                                  &tcpsession::handler
                                , this
                                , boost::asio::placeholders::error
                                , boost::asio::placeholders::bytes_transferred
                            )
                        );
                    else
                        delete this;
                }

            private:

                void handler(const boost::system::error_code& error, boost::uint32_t size) {
                    if (!error) {
                        std::string action(data, size);
                        std::string response;
                        if (actions.find(action) != actions.end())
                            response = actions[action]();
                        else
                            response = "Unknown action: " + action;
                        std::memcpy(data + 4, response.c_str(), size = response.size());
                        std::memcpy(data, &size, 4);
                        boost::asio::async_write(
                              socket
                            , boost::asio::buffer(data, 1028)
                            , boost::bind(&tcpsession::start, this, boost::asio::placeholders::error)
                        );
                    } else
                        delete this;
                }

                actions_type & actions;
                char data[1028];
                boost::asio::ip::tcp::socket socket;
        };
    }
}

#endif
