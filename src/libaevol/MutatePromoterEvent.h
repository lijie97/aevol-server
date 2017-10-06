//
// Created by Jonathan Rouzaud-Cornabas on 24/08/17.
//

#ifndef RAEVOL_CUDA_MUTATEPROMOTEREVENT_H
#define RAEVOL_CUDA_MUTATEPROMOTEREVENT_H

#include <cstdint>
#include "SIMD_Individual.h"

namespace aevol {


enum MutatePromoterEventType {
    REMOVE           = 0,
    REMOVE_LEADING      = 1,
    REMOVE_LAGGING = 2,
    LOCATE               = 3,
    LOCATE_LEADING         = 4,
    LOCATE_LAGGING        = 5,
    MOVE_LEADING         = 6,
    MOVE_LAGGING         = 7,
    DUPLICATE            = 8,
    DO_INSERT               = 9,
    EXTRACT_LEADING      = 10,
    EXTRACT_LAGGING      = 11,
    INVERT               = 12,
    INSERT_POS           = 13,
    SHIFT_POS            = 14
};


class MutatePromoterEvent {
 public:
    MutatePromoterEvent() = default;

    void remove_all();
    void remove_leading(int32_t pos1, int32_t pos2);
    void remove_lagging(int32_t pos1, int32_t pos2);

    void locate_all();
    void locate_leading(int32_t pos1, int32_t pos2);
    void locate_lagging(int32_t pos1, int32_t pos2);

    void move_leading(int32_t pos1, int32_t pos2);
    void move_lagging(int32_t pos1, int32_t pos2);

    void duplicate(int32_t pos1, int32_t pos2);

    void insert(int32_t pos);
    void insert();

    void extract_leading(int32_t pos1, int32_t pos2);
    void extract_lagging(int32_t pos1, int32_t pos2);

    void shift_pos(int32_t pos1, int32_t pos2);

    void invert(int32_t pos1, int32_t pos2);

    int32_t type_;
    bool all_;
    int32_t pos_1_,pos_2_;

    //list<promoterStruct> tmp_list;
};
}


#endif //RAEVOL_CUDA_MUTATEPROMOTEREVENT_H
