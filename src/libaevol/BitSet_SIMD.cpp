//
// Created by arrouan on 24/01/18.
//

#include <mm_malloc.h>
#include <cstring>
#include <cstdio>
#include "BitSet_SIMD.h"

aevol::BitSet_SIMD::BitSet_SIMD(const char* data, int32_t length) {
  length_ = length;

  nb_char_block_ = (length_/ sizeof(char)) + 1;

  nb_block_ = nb_char_block_/BITSET_BLOCK_SIZE + 1;

  posix_memalign((void **)&data_ ,64, nb_block_ * BITSET_BLOCK_SIZE* sizeof(char));

  memset(data_,0,nb_char_block_);

  for (int32_t i = 0; i < length_; i++) {
    if (data[i]=='0')
      set_to_0(i);
    else
      set_to_1(i);
  }
}

aevol::BitSet_SIMD::BitSet_SIMD(aevol::BitSet_SIMD* bitset) {
  length_ = bitset->length_;
  nb_char_block_ = bitset->nb_char_block_;
  nb_block_ = bitset->nb_block_;

  posix_memalign((void **)&data_ ,64, nb_block_ * BITSET_BLOCK_SIZE* sizeof(char));

  memcpy(data_, bitset->data_, (nb_char_block_) * sizeof(char));
}

void aevol::BitSet_SIMD::print() {
  printf("Length %d\n",length_);
  for (int32_t i = 0; i < length_; i++) {
    printf("\n %d ",i);
    if (get(i))
      printf("1");
    else
      printf("0");
  }
  printf("\n");
}

void aevol::BitSet_SIMD::append(aevol::BitSet_SIMD* bitset) {
  char* new_data;
  int old_length = length_;

  if (length_+bitset->length_ > nb_block_*BITSET_BLOCK_SIZE*sizeof(char)) {
    length_ = length_+bitset->length_;

    nb_char_block_ = (length_/ sizeof(char)) + 1;

    nb_block_ = nb_char_block_/BITSET_BLOCK_SIZE + 1;

    posix_memalign((void **)&new_data ,64, nb_block_ * BITSET_BLOCK_SIZE* sizeof(char));
    memcpy(new_data, data_, (nb_char_block_) * sizeof(char));
    free(data_);
    data_ = new_data;
  } else
    length_ = length_+bitset->length_;

  for (int i = old_length; i < length_; i++)
    bitset->get(i-old_length) ? set_to_1(i) : set_to_0(i);
}

void aevol::BitSet_SIMD::insert_at(aevol::BitSet_SIMD* bitset, int position) {
  char* new_data;
  int old_length = length_;
  char* old_data = data_;

  if (length_+bitset->length_ > nb_block_*BITSET_BLOCK_SIZE*sizeof(char)) {
    length_ = length_+bitset->length_;

    nb_char_block_ = (length_/ sizeof(char)) + 1;

    nb_block_ = nb_char_block_/BITSET_BLOCK_SIZE + 1;

    posix_memalign((void **)&new_data ,64, nb_block_ * BITSET_BLOCK_SIZE* sizeof(char));
    memcpy(new_data, data_, (nb_char_block_) * sizeof(char));
    data_ = new_data;

  } else
    length_ = length_+bitset->length_;

  for (int i = position+bitset->length_; i < length_; i++) {
    old_data[((i - bitset->length_) / sizeof(char))] &
    (1 << ((i - bitset->length_) % sizeof(char))) ?
    set_to_1(i) : set_to_0(i);
  }

  for (int i = position; i < position+bitset->length_; i++)
    bitset->get(i-position) ? set_to_1(i) : set_to_0(i);


  free(old_data);
}

/*
int main() {
  const char data[] = "110111";

  printf("CHAR[] : %s\n",data);
  aevol::BitSet_SIMD* bitset = new aevol::BitSet_SIMD(data,6);
  bitset->print();

  aevol::BitSet_SIMD* clone_bitset = new aevol::BitSet_SIMD(bitset);
  clone_bitset->print();

  bitset->flip(1);
  bitset->print();

  //bitset->append(clone_bitset);
  bitset->insert_at(clone_bitset,3);
  bitset->print();
}*/
