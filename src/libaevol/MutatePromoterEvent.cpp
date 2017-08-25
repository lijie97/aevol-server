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

void MutatePromoterEvent::duplicate(int32_t pos1, int32_t pos2,
                                    list<promoterStruct> list) {
  type_ = MutatePromoterEventType::DUPLICATE;
  pos_1_ = pos1;
  pos_2_ = pos2;
  tmp_list = list;
}

void MutatePromoterEvent::insert(std::list<promoterStruct> list,  int32_t pos) {
  type_ = MutatePromoterEventType::INSERT_POS;
  pos_1_ = pos;
  tmp_list = list;
}

void MutatePromoterEvent::insert(std::list<promoterStruct> list) {
  type_ = MutatePromoterEventType::INSERT;
  tmp_list = list;
}

void MutatePromoterEvent::extract_leading(int32_t pos1, int32_t pos2,
                                          std::list <promoterStruct> list) {
  type_ = MutatePromoterEventType::EXTRACT_LEADING;
  pos_1_ = pos1;
  pos_2_ = pos2;
  tmp_list = list;
}


void MutatePromoterEvent::extract_lagging(int32_t pos1, int32_t pos2,
                                          std::list <promoterStruct> list) {
  type_ = MutatePromoterEventType::EXTRACT_LAGGING;
  pos_1_ = pos1;
  pos_2_ = pos2;
  tmp_list = list;
}

void MutatePromoterEvent::invert(std::list<promoterStruct> list, int32_t pos1,
                                 int32_t pos2) {
  type_ = MutatePromoterEventType::INVERT;
  pos_1_ = pos1;
  pos_2_ = pos2;
}
}
