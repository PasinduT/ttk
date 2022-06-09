#include <Testing.h>

ttk::Testing::Testing() {
  // inherited from Debug: prefix will be printed at the beginning of every msg
  this->setDebugMsgPrefix("Testing");
}
