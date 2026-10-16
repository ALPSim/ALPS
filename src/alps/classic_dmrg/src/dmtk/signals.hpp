// Ported unchanged from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file
// src/solvers/tensor_network/src/classic_dmrg/dmtk/signals.hpp, itself identical to
// src/apps/dmrg/dmrg/dmtk/signals.h on master (Adrian Feiguin's dmtk), which
// carries no notice of its own.

#ifndef __DMTK_SIGNALS__
#define __DMTK_SIGNALS__

namespace dmtk
{

template<class A>
class Signal
{
  public:
    typedef bool (*callback) (A &object, size_t signal_id, void *data);

    Signal() {}
    Signal(callback handler, A *object, size_t signal_id, void *data):
      _handler(handler), _object(object), _signal_id(signal_id), _data(data) {}

    void * data() { return _data; }
    const void * data() const { return _data; }
    size_t & signal_id() { return _signal_id; }
    size_t signal_id() const { return _signal_id; }

    bool emit()
      { return _handler(*_object, _signal_id, _data); }

  private:
    callback _handler;
    A *_object;
    size_t _signal_id;
    void *_data;
};


} // namespace dmtk

#endif // __DMTK_SIGNALS__
