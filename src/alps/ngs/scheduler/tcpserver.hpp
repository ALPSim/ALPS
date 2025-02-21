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

#ifndef ALPS_NGS_TCPSERVER_HPP
#define ALPS_NGS_TCPSERVER_HPP

#include <alps/ngs/detail/tcpsession.hpp>

#include <boost/bind.hpp>
#include <boost/asio.hpp>
#include <boost/function.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/lexical_cast.hpp>

#include <map>
#include <string>

namespace alps {

    class tcpserver {
        
        public:
        
            typedef std::map<std::string, boost::function<std::string()> > actions_type;

            tcpserver(short port, actions_type const & a = actions_type())
                : actions(a)
                , io_service()
                , acceptor(io_service, boost::asio::ip::tcp::endpoint(boost::asio::ip::tcp::v4(), port))
            {
                start();
            }
        
            void add_action(std::string const & name, boost::function<std::string()> const & action) {
                actions[name] = action;
            }
        
            void poll() {
                io_service.poll();
            }
        
            void stop() {
                io_service.stop();
            }

        private:

            void start() {
                detail::tcpsession * new_sess = new detail::tcpsession(actions, io_service); // TODO: use shared_ptr and manage the sockets manually
                acceptor.async_accept(
                      new_sess->get_socket()
                    , boost::bind(&tcpserver::handler, this, new_sess, boost::asio::placeholders::error)
                );
            }

            void handler(detail::tcpsession * new_sess, boost::system::error_code const & error) {
                if (!error)
                    new_sess->start();
                else
                    delete new_sess;
                start();
            }

            actions_type actions;
            boost::asio::io_service io_service;
            boost::asio::ip::tcp::acceptor acceptor;
    };

}

#endif
