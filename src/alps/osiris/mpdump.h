/***************************************************************************
* PALM++/osiris library
*
* osiris/mpdump.h   message passing via dumps
*
* $Id$
*
* Copyright (C) 1994-2002 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
*                            Synge Todo <wistaria@comp-phys.org>,
*
* This program is free software; you can redistribute it and/or
* modify it under the terms of the GNU General Public License
* as published by the Free Software Foundation; either version 2
* of the License, or (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
**************************************************************************/

//=======================================================================
// This file contains the classes for message passsing
// messages are objects dumped int32_to a message passing dump
// and sent to a process, where they can be reconstructed from the dump
//=======================================================================

#ifndef OSIRIS_MPDUMP_H
#define OSIRIS_MPDUMP_H

#include <alps/config.h>
#include <alps/osiris/dump.h>
#include <alps/osiris/process.h>
#include <alps/osiris/buffer.h>

namespace alps {

//=======================================================================
// alps::OMPDump
//
// objects written into this dump can be sent to another process, where
// they can be received as an alps::IDump
//-----------------------------------------------------------------------

class OMPDump : public ODump
{
public:
  OMPDump(); 
  ~OMPDump(); 

  // MESSAGE PASSING FUNCTIONS
  
  // initialize the send buffer
  void init(); 

  void send(const Process&,int32_t tag); // send to a single process
  void send(const ProcessList&,int32_t); // broadcast to a whole group

# define ALPS_DUMP_DO_TYPE(T) \
  void write_simple(T x); \
  void write_array(std::size_t, const T *);
  ALPS_DUMP_DO_TYPE(bool)
  ALPS_DUMP_DO_TYPE(char)
  ALPS_DUMP_DO_TYPE(signed char)
  ALPS_DUMP_DO_TYPE(unsigned char)
  ALPS_DUMP_DO_TYPE(short)
  ALPS_DUMP_DO_TYPE(unsigned short)
  ALPS_DUMP_DO_TYPE(int)
  ALPS_DUMP_DO_TYPE(unsigned int)
  ALPS_DUMP_DO_TYPE(long)
  ALPS_DUMP_DO_TYPE(unsigned long)
# ifdef BOOST_HAS_LONG_LONG
  ALPS_DUMP_DO_TYPE(long long)
  ALPS_DUMP_DO_TYPE(unsigned long long)
# endif
  ALPS_DUMP_DO_TYPE(float)
  ALPS_DUMP_DO_TYPE(double)
  ALPS_DUMP_DO_TYPE(long double)
# undef ALPS_DUMP_DO_TYPE
    
  // write a c-style string
  void write_string(std::size_t, const char *);

private:
  bool valid_; // flag to indicate the state of the message buffers
  int bufid_; // the PVM buffer id
  int oldbufid_; //the previously used buffer
  detail::Buffer buf_; // the message buffer for other message passing systems
};


//=======================================================================
// IMPDump
//
// derived from alps::IDump
//
// can be received as a message
//-----------------------------------------------------------------------

class IMPDump : public IDump
{
public:
  
   static int32_t probe(int32_t=-1);
   static int32_t probe(const Process& w, int32_t t=-1); 

   IMPDump(); 
   IMPDump(int32_t t); 
   IMPDump(const Process&, int32_t t); 

  // reinitialize the buffer
  void init();

  const Process& sender() const;   // the sender of the message
  void receive(const Process& w,int32_t t); // receive a message from a specific process
  void receive(int32_t t); // receive a mesage from anywhere

# define ALPS_DUMP_DO_TYPE(T) \
  void read_simple(T& x); \
  void read_array(std::size_t, T *);
  ALPS_DUMP_DO_TYPE(bool)
  ALPS_DUMP_DO_TYPE(char)
  ALPS_DUMP_DO_TYPE(signed char)
  ALPS_DUMP_DO_TYPE(unsigned char)
  ALPS_DUMP_DO_TYPE(short)
  ALPS_DUMP_DO_TYPE(unsigned short)
  ALPS_DUMP_DO_TYPE(int)
  ALPS_DUMP_DO_TYPE(unsigned int)
  ALPS_DUMP_DO_TYPE(long)
  ALPS_DUMP_DO_TYPE(unsigned long)
# ifdef BOOST_HAS_LONG_LONG
  ALPS_DUMP_DO_TYPE(long long)
  ALPS_DUMP_DO_TYPE(unsigned long long)
# endif
  ALPS_DUMP_DO_TYPE(float)
  ALPS_DUMP_DO_TYPE(double)
  ALPS_DUMP_DO_TYPE(long double)
# undef ALPS_DUMP_DO_TYPE
    
  // read a c-style string
  void read_string(std::size_t, char *);

private:
  bool valid_; // state of the message buffer
  int bufid_; // the PVM buffer id
  detail::Buffer buf_; // the message buffer
  Process theSender_; // the process that sent the message
  void receive(const Process*,int32_t);   // receive a message
  static int32_t probe(const Process* w, int32_t t=-1); 
};

}

#endif // OSIRIS_MPDUMP_H
