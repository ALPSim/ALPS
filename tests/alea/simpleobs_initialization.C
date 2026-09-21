// SPDX-License-Identifier: MIT

#include <alps/alea.h>
#include <cstring>
#include <iostream>
#include <new>

int main() {
  using data_type = alps::SimpleObservableData<double>;
  alignas(data_type) unsigned char storage[sizeof(data_type)];
  // Reused storage must not change the default observable's public state.
  std::memset(storage, 1, sizeof(storage));
  auto* data = new (storage) data_type;
  const bool can_thermalize = data->can_set_thermalization();
  data_type copy(*data);
  data->~data_type();
  if (can_thermalize || copy.can_set_thermalization()) {
    std::cerr << "A default observable must not allow setting thermalization\n";
    return 1;
  }
}
