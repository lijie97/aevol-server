//
// Created by Jonathan Rouzaud-Cornabas on 24/08/17.
//

#include "MutationEvent.h"

namespace aevol {

void MutationEvent::switch_pos(int32_t pos) {
  type_ = MutationEventType::DO_SWITCH;
  pos_1_ = pos;
}

#ifdef WITH_BITSET
void MutationEvent::small_insertion(int32_t pos, BitSet_SIMD* seq) {
  type_ = MutationEventType::SMALL_INSERTION;
  pos_1_ = pos;
  seq_ = seq;
}
#else
void MutationEvent::small_insertion(int32_t pos, int32_t number, char* seq) {
  type_ = MutationEventType::SMALL_INSERTION;
  pos_1_ = pos;
  number_ = number;
  seq_ = seq;
}
#endif

void MutationEvent::small_deletion(int32_t pos, int32_t number) {
  type_ = MutationEventType::SMALL_DELETION;
  pos_1_ = pos;
  number_ = number;
}

void MutationEvent::duplication(int32_t pos1, int32_t pos2, int32_t pos3) {
  type_ = MutationEventType::DUPLICATION;
  pos_1_ = pos1;
  pos_2_ = pos2;
  pos_3_ = pos3;
}

void MutationEvent::translocation(int32_t pos1, int32_t pos2, int32_t pos3,
                                    int32_t pos4, bool invert) {
  type_ = MutationEventType::TRANSLOCATION;
  pos_1_ = pos1;
  pos_2_ = pos2;
  pos_3_ = pos3;
  pos_4_ = pos4;
  invert_ = invert;
}

void MutationEvent::inversion(int32_t pos1, int32_t pos2) {
  type_ = MutationEventType::INVERSION;
  pos_1_ = pos1;
  pos_2_ = pos2;
}

void MutationEvent::deletion(int32_t pos1, int32_t pos2) {
  type_ = MutationEventType::DELETION;
  pos_1_ = pos1;
  pos_2_ = pos2;
}

MutationEvent::~MutationEvent() {
  if (seq_ != nullptr)
#ifdef WITH_BITSET
    delete seq_;
#else
    delete [] seq_;
#endif
}
}
