//
// Created by Jonathan Rouzaud-Cornabas on 24/08/17.
//

#ifndef RAEVOL_CUDA_MUTATIONEVENT_H
#define RAEVOL_CUDA_MUTATIONEVENT_H


#include <cstdint>

namespace aevol {

enum MutationEventType {
    DO_SWITCH           = 0,
    SMALL_INSERTION      = 1,
    SMALL_DELETION = 2,
    DUPLICATION               = 3,
    TRANSLOCATION         = 4,
    INVERSION        = 5,
    DELETION         = 6
};

class MutationEvent {

 public:
    MutationEvent() = default;
    ~MutationEvent();

    void switch_pos(int32_t pos);
    void small_insertion(int32_t pos, int32_t number, char* seq);
    void small_deletion(int32_t pos, int32_t number);

    void duplication(int32_t pos1, int32_t pos2, int32_t pos3);

    void translocation(int32_t pos1, int32_t pos2, int32_t pos3,
                       int32_t pos4, bool invert);

    void inversion(int32_t pos1, int32_t pos2);

    void deletion(int32_t pos1, int32_t pos2);

    int32_t type() { return type_; };

    int32_t pos_1() { return pos_1_; }
    int32_t pos_2() { return pos_2_; }
    int32_t pos_3() { return pos_3_; }
    int32_t pos_4() { return pos_4_; }

    int32_t number() { return number_; }

    int32_t invert() { return invert_; }

    char* seq() { return seq_; }

 private:
    int32_t type_;

    int32_t pos_1_,pos_2_,pos_3_,pos_4_;

    int32_t number_; // insertion or deletion

    bool invert_;

    char* seq_ = nullptr;
};
}


#endif //RAEVOL_CUDA_MUTATIONEVENT_H
