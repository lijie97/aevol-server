//
// Created by Jonathan Rouzaud-Cornabas on 24/08/17.
//

#include "MutatePromoterEvent.h"

namespace aevol {
void MutatePromoterEvent::remove_all() {
  type_ = MutatePromoterEventType::REMOVE;
  all_ = true;
}

void MutatePromoterEvent::remove_leading(int32_t pos1, int32_t pos2) {
  type_ = MutatePromoterEventType::REMOVE_LEADING;
  pos_1_ = pos1;
  pos_2_ = pos2;
}

void MutatePromoterEvent::remove_lagging(int32_t pos1, int32_t pos2) {
  type_ = MutatePromoterEventType::REMOVE_LAGGING;
  pos_1_ = pos1;
  pos_2_ = pos2;
}

void MutatePromoterEvent::locate_all() {
  type_ = MutatePromoterEventType::LOCATE;
  all_ = true;
}

void MutatePromoterEvent::locate_leading(int32_t pos1, int32_t pos2) {
  type_ = MutatePromoterEventType::LOCATE_LEADING;
  pos_1_ = pos1;
  pos_2_ = pos2;
}

void MutatePromoterEvent::locate_lagging(int32_t pos1, int32_t pos2) {
  type_ = MutatePromoterEventType::LOCATE_LAGGING;
  pos_1_ = pos1;
  pos_2_ = pos2;
}

void MutatePromoterEvent::move_leading(int32_t pos1, int32_t pos2) {
  type_ = MutatePromoterEventType::MOVE_LEADING;
  pos_1_ = pos1;
  pos_2_ = pos2;
}

void MutatePromoterEvent::move_lagging(int32_t pos1, int32_t pos2) {
  type_ = MutatePromoterEventType::MOVE_LAGGING;
  pos_1_ = pos1;
  pos_2_ = pos2;
}

void MutatePromoterEvent::duplicate(int32_t pos1, int32_t pos2) {
  type_ = MutatePromoterEventType::DUPLICATE;
  pos_1_ = pos1;
  pos_2_ = pos2;
}

void MutatePromoterEvent::insert(int32_t pos) {
  type_ = MutatePromoterEventType::INSERT_POS;
  pos_1_ = pos;
}

void MutatePromoterEvent::insert() {
  type_ = MutatePromoterEventType::DO_INSERT;
}

void MutatePromoterEvent::extract_leading(int32_t pos1, int32_t pos2) {
  type_ = MutatePromoterEventType::EXTRACT_LEADING;
  pos_1_ = pos1;
  pos_2_ = pos2;
}


void MutatePromoterEvent::extract_lagging(int32_t pos1, int32_t pos2) {
  type_ = MutatePromoterEventType::EXTRACT_LAGGING;
  pos_1_ = pos1;
  pos_2_ = pos2;
}

void MutatePromoterEvent::invert(int32_t pos1,
                                 int32_t pos2) {
  type_ = MutatePromoterEventType::INVERT;
  pos_1_ = pos1;
  pos_2_ = pos2;
}

void MutatePromoterEvent::shift_pos(int32_t pos1, int32_t pos2) {
  type_ = MutatePromoterEventType::SHIFT_POS;
  pos_1_ = pos1;
  pos_2_ = pos2;
}
}
