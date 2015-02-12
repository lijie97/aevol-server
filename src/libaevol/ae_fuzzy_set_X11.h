// Aevol - An in silico experimental evolution platform
// 
// Copyright: See the AUTHORS file provided with the package or <www.aevol.fr>
// Web: http://www.aevol.fr/
// E-mail: See <http://www.aevol.fr/contact/>
// Original Authors : Guillaume Beslon, Carole Knibbe, David Parsons
// 
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
 
#ifndef __AE_FUZZY_SET_X11_H__
#define __AE_FUZZY_SET_X11_H__
 
#include <inttypes.h>
#include <stdlib.h>
#include <stdio.h>

#include "fuzzy.h"
#include <ae_X11_window.h>

namespace aevol {
class ae_fuzzy_set_X11 : public Fuzzy
{  
 public :
  ae_fuzzy_set_X11(): Fuzzy() {};
  ae_fuzzy_set_X11(const Fuzzy& model): Fuzzy(model) {};
  virtual ~ae_fuzzy_set_X11() {};  
  virtual void display(ae_X11_window* win, color_map color, bool fill = false, bool bold = false);
};
} // namespace aevol
#endif // __AE_FUZZY_SET_X11_H__
