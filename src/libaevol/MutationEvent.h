//
// Created by Jonathan Rouzaud-Cornabas on 24/08/17.
//

#ifndef RAEVOL_CUDA_MUTATIONEVENT_H
#define RAEVOL_CUDA_MUTATIONEVENT_H


#include <cstdint>
#include <set>
#include "MutatePromoterEvent.h"

namespace aevol {

enum MutationEventType {
    SWITCH           = 0,
    SMALL_INSERTION      = 1,
    SMALL_DELETION = 2,
    DUPLICATION               = 3,
    TRANSLOCATION         = 4,
    INVERSION        = 5,
    DELETION         = 6
};

class MutationEvent {

 public:
    MutationEvent() = default {};

    void switch_pos(int32_t pos);
    void small_insertion(int32_t pos, int32_t number);
    void small_deletion(int32_t pos, int32_t number);

    void duplication(int32_t pos1, int32_t pos2, int32_t pos3);

    void translocation(int32_t pos1, int32_t pos2, int32_t pos3,
                       int32_t pos4, bool invert);

    void inversion(int32_t pos1, int32_t pos2);

    void deletion(int32_t pos1, int32_t pos2);

 private:
    int32_t type_;

    int32_t pos_1_,pos_2_,pos_3_,pos_4_;

    int32_t number_; // insertion or deletion

    bool invert_;

    std::set<MutatePromoterEvent> promoter_event_list_;
};
}


#endif //RAEVOL_CUDA_MUTATIONEVENT_H
