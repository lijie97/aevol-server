//
// Created by arrouan on 24/01/18.
//

#include <mm_malloc.h>
#include <cstring>
#include <cstdio>
#include <cassert>
#include "BitSet_SIMD.h"

aevol::BitSet_SIMD::BitSet_SIMD(const char* data, int32_t length) {
  length_ = length;

#ifdef _DYNAMIC_CUSTOM_BITSET
  nb_char_block_ = (length_/ (8*sizeof(char))) + 1;

  nb_block_ = nb_char_block_/BITSET_BLOCK_SIZE + 1;

  //printf("Length %d Char block %d block %d\n",length_,nb_char_block_,nb_block_);

  //printf("%d\n",nb_block_ * BITSET_BLOCK_SIZE* sizeof(char));

  //data_ = (char *) calloc(nb_block_ * BITSET_BLOCK_SIZE ,sizeof(char));
  posix_memalign((void **)&data_ ,32, nb_block_ * BITSET_BLOCK_SIZE* sizeof(char));
  //if (data == NULL) printf("NNNNN\n\n");

  memset(data_,0,nb_block_ * BITSET_BLOCK_SIZE* sizeof(char));

  for (int32_t i = 0; i < length_; i++) {
    //printf("Loading %d\n",i);
    if (data[i]=='0')
      set_to_0(i);
    else
      set_to_1(i);
  }
#elif _STATIC_BITSET
  for (int32_t i = 0; i < length_; i++) {
    //printf("Loading %d\n",i);
    if (data[i]=='1')
      data_.set(i);
  }
#endif

  //printf("Going out\n");
}

aevol::BitSet_SIMD::BitSet_SIMD(aevol::BitSet_SIMD* bitset) {
  length_ = bitset->length_;
#ifdef _DYNAMIC_CUSTOM_BITSET
  nb_char_block_ = bitset->nb_char_block_;
  nb_block_ = bitset->nb_block_;

  posix_memalign((void **)&data_ ,64, nb_block_ * BITSET_BLOCK_SIZE* sizeof(char));

  memcpy(data_, bitset->data_, (nb_char_block_) * sizeof(char));
#elif _STATIC_BITSET
  data_ = bitset->data_;
#endif
}

aevol::BitSet_SIMD::BitSet_SIMD(int32_t length) {
  length_ = length;
#if _DYNAMIC_CUSTOM_BITSET
  nb_char_block_ = (length_/ (8*sizeof(char))) + 1;

  nb_block_ = nb_char_block_/BITSET_BLOCK_SIZE + 1;

  posix_memalign((void **)&data_ ,64, nb_block_ * BITSET_BLOCK_SIZE* sizeof(char));
  memset(data_,0,nb_block_ * BITSET_BLOCK_SIZE* sizeof(char));
#elif _STATIC_BITSET

#endif
}

void aevol::BitSet_SIMD::print() {
  printf("Length %d\n",length_);
  for (int32_t i = 0; i < length_; i++) {
    printf("\n %d ",i);
#ifdef _DYNAMIC_CUSTOM_BITSET
    if (get(i))
#elif _STATIC_BITSET
    if (data_[i])
#else
    if (true)
#endif
      printf("1");
    else
      printf("0");
  }
  printf("\n");
}

void aevol::BitSet_SIMD::append(aevol::BitSet_SIMD* bitset) {
  int32_t old_length = length_;
#ifdef _DYNAMIC_STATIC_BITSET
  char* new_data;


  if (length_+bitset->length_ > nb_block_*BITSET_BLOCK_SIZE*sizeof(char)) {
    length_ = length_+bitset->length_;

    int old_char_block = nb_char_block_;
    nb_char_block_ = (length_/ (8*sizeof(char))) + 1;

    nb_block_ = nb_char_block_/BITSET_BLOCK_SIZE + 1;

    posix_memalign((void **)&new_data ,64, nb_block_ * BITSET_BLOCK_SIZE* sizeof(char));
    memcpy(new_data, data_, (old_char_block) * sizeof(char));
    free(data_);
    data_ = new_data;
  } else
    length_ = length_+bitset->length_;

  for (int32_t i = old_length; i < length_; i++)
    bitset->get(i-old_length) ? set_to_1(i) : set_to_0(i);
#elif _STATIC_BITSET
  length_ = length_+bitset->length_;

  for (int32_t i = old_length; i < length_; i++)
    data_.set(i,bitset->data_[i-old_length]);
#endif
}

void aevol::BitSet_SIMD::insert_at(aevol::BitSet_SIMD* bitset, int32_t position) {

  int32_t old_length = length_;


#ifdef _DYNAMIC_CUSTOM_BITSET
  char* new_data;

  char* old_data = data_;
  bool to_free = false;

  if (length_+bitset->length_ > nb_block_*BITSET_BLOCK_SIZE*sizeof(char)) {
    length_ = length_+bitset->length_;

    int old_char_block = nb_char_block_;
    nb_char_block_ = (length_/ (8*sizeof(char))) + 1;

    nb_block_ = nb_char_block_/BITSET_BLOCK_SIZE + 1;

    posix_memalign((void **)&new_data ,64, nb_block_ * BITSET_BLOCK_SIZE* sizeof(char));

    memcpy(new_data, data_, (old_char_block) * sizeof(char));
    data_ = new_data;
    to_free = true;
  } else {
    length_ = length_ + bitset->length_;
  }

  for (int32_t i = position+bitset->length_; i < length_; i++) {
    get(old_data,i - bitset->length_) ?
    set_to_1(i) : set_to_0(i);
  }

  for (int32_t i = position; i < position+bitset->length_; i++)
    bitset->get(i-position) ? set_to_1(i) : set_to_0(i);


  if (to_free)
    free(old_data);
#elif _STATIC_BITSET
  length_ = length_ + bitset->length_;

  std::bitset<BITSET_STATIC_MAX_SIZE> tmp = data_;

  for (int32_t i = position+bitset->length_; i < length_; i++) {
    tmp.set(i,data_[i - bitset->length_]);
  }

  for (int32_t i = position; i < position+bitset->length_; i++)
    tmp.set(i-position,bitset->data_[i-position]);

  data_ = tmp;
#endif
}

void aevol::BitSet_SIMD::remove(int32_t pos_1, int32_t pos_2) {
  assert(pos_1 >= 0 && pos_2 >= pos_1 && pos_2 <= length_);


  int32_t old_length = length_;
  length_    = length_ - (pos_2 - pos_1);

#ifdef _DYNAMIC_CUSTOM_BITSET
  char* new_data;
  char* old_data = data_;


  // Compute size of new genome
  int old_char_block = nb_char_block_;

  nb_char_block_ = (length_/ (8*sizeof(char))) + 1;

  nb_block_ = nb_char_block_/BITSET_BLOCK_SIZE + 1;

  posix_memalign((void **)&new_data ,64, nb_block_ * BITSET_BLOCK_SIZE* sizeof(char));
  memset(new_data,0,nb_block_ * BITSET_BLOCK_SIZE* sizeof(char));

  // Copy the remaining of the genome in tmp (preceeding and following parts)
  memcpy(new_data, data_, ((int)(pos_1/8)+1) * sizeof(char));

  data_ = new_data;

  //printf("Reduce from %d to %d :: pos1 %d pos2 %d (copy from 0 to %d)\n",old_length,length_,pos_1,pos_2,old_length-pos_2);

  for (int32_t i = 0; i < old_length-pos_2; i++) {
    get(old_data,pos_2+i) ?
    set_to_1(pos_1+i) : set_to_0(pos_1+i);
  }

  free(old_data);
#elif _STATIC_BITSET
  std::bitset<BITSET_STATIC_MAX_SIZE> tmp;

  tmp = data_;

  for (int32_t i = 0; i < old_length-pos_2; i++) {
    tmp.set(pos_1+i,data_[pos_2+i]);
  }
#endif
}

aevol::BitSet_SIMD* aevol::BitSet_SIMD::duplicate(int32_t pos_1, int32_t pos_2, bool invert) {
  BitSet_SIMD* duplicate_bitset = new BitSet_SIMD(pos_2-pos_1);

  if (invert) {
    for (int32_t i = 0, j = pos_2 - 1; i < duplicate_bitset->length_; i++, j--) {
      //printf("%d %d\n",i,j);
#ifdef _DYNAMIC_CUSTOM_BITSET
      get(j) ? duplicate_bitset->set_to_1(i) : duplicate_bitset->set_to_0(i);
#elif _STATIC_BITSET
      duplicate_bitset->data_.set(i,data_[j]);
#endif
    }
  } else {
    for (int32_t i = 0; i < duplicate_bitset->length_; i++) {
#ifdef _DYNAMIC_CUSTOM_BITSET
      get(pos_1+i) ?
      duplicate_bitset->set_to_1(i) : duplicate_bitset->set_to_0(i);
#elif _STATIC_BITSET
      duplicate_bitset->data_.set(i,data_[pos_1+i]);
#endif
    }
  }

  return duplicate_bitset;
}

void aevol::BitSet_SIMD::replace(aevol::BitSet_SIMD* bitset, int32_t position) {
  //printf("%d to %d by %d\n",position,length_,bitset->length_);


  for (int32_t i = position; i < position+bitset->length_; i++) {
#ifdef _DYNAMIC_CUSTOM_BITSET
    bitset->get(i-position) ? set_to_1(i) : set_to_0(i);
#elif _STATIC_BITSET
    data_.set(i,bitset->data_[i-position]);
#endif
  }
}


void aevol::BitSet_SIMD::replace(aevol::BitSet_SIMD* bitset, int32_t start, int32_t stop, int32_t position) {
  for (int32_t i = position, j = start; i < position+(stop-start); i++, j++) {
#ifdef _DYNAMIC_CUSTOM_BITSET
    bitset->get(j) ? set_to_1(i) : set_to_0(i);
#elif _STATIC_BITSET
    data_.set(i,bitset->data_[j]);
#endif
  }
}

void aevol::BitSet_SIMD::insert_at(aevol::BitSet_SIMD* bitset, int32_t start,
                                   int32_t stop, int32_t position) {

  int32_t old_length = length_;
  int32_t to_add_length = stop - start;

#ifdef _DYNAMIC_CUSTOM_BITSET
  char* old_data = data_;
  char* new_data;


  bool to_free = false;
  //printf("insert at %d : start %d stop %d length %d current length %d\n",position,start,stop,length_+to_add_length,length_);

  if (length_+to_add_length > nb_block_*BITSET_BLOCK_SIZE*sizeof(char)) {
    length_ = length_+to_add_length;

    int old_char_block = nb_char_block_;
    nb_char_block_ = (length_/ (8*sizeof(char))) + 1;

    nb_block_ = nb_char_block_/BITSET_BLOCK_SIZE + 1;

    posix_memalign((void **)&new_data ,64, nb_block_ * BITSET_BLOCK_SIZE* sizeof(char));

    memset(new_data,0,nb_block_ * BITSET_BLOCK_SIZE* sizeof(char));

    memcpy(new_data, data_, (old_char_block) * sizeof(char));
    data_ = new_data;
    to_free = true;
    //printf("New size is %d char block %d nb block %d\n",length_,nb_char_block_,nb_block_);
  } else {
    length_ = length_+to_add_length;
  }

  for (int32_t i = position; i < old_length; i++) {
    //printf("MOVING %d to %d\n",i,i+to_add_length);
    get(old_data,i) ?
      set_to_1(i+to_add_length) : set_to_0(i+to_add_length);
  }

  for (int32_t i = start; i < stop; i++)
    bitset->get(i) ? set_to_1(i+position) : set_to_0(i+position);


  if (to_free)
    free(old_data);
#elif _STATIC_BITSET
  length_ = length_+to_add_length;

  std::bitset<BITSET_STATIC_MAX_SIZE> tmp = data_;

  for (int32_t i = position; i < old_length; i++) {
    tmp.set(i+to_add_length,bitset->data_[i]);
  }

  for (int32_t i = start; i < stop; i++)
    tmp.set(i+position,bitset->data_[i]);

  data_ = tmp;

  delete bitset;
#endif
}

/*
int main() {
  const char data[] = "01001101111110000101011001110010010110011010010011100110101001010100101100010001111010111110001101011000111111111000010100110110010010001100000011011001101011101100000010101011001101100100001011000101100011101100111000100000110000101110100010110101111110001011110011110110101010110101111100110101111001101101111110111111111100100001111000001101100111111100111101001110101100010000110000100011000110110010001011110000000111110000011010001101000111000011110110010011100110111000101110000111000010110000111100110101001100011101111111001111010001100100010101011111100001000001001001001000010011010111101111111111100101111101100011000100011100010001000100101101110010010110100111111010000101000001000010000011110001000010000111001011100001000111100011011111010010110111001001100001000101110001111011101010101011110111111010110000010010000001111101101001111001000111111111101110110000000000110111000011111001101010101011111101011011100100011011110011000011010100100111001101011000111111010111000101010101110101000100110010100100111001100001001100100010011100001001001110100110011110110011011100100101011000100110100101000110011110101010101011111101001101000101000001011110100000011000111010001011001101000010111110111000011110001110010101001111111100111001100011100111110111011011111000000000110111011001101110110100101011010110100000011001100101110000110010111001000001111010100111001101010111000011110000110001010100000000000100100100111000010010110101110111011101111001101000001100111100000010010010111001100110010100011111001011010101110111001101010001011111101000100000101010000000101010100011011011000100100101100101000001110100011100000011100100010101011101111010010111101111110000000000000011011000000001000110111101001000111101101110111000101111001011111100011101101111101101001101111010001000101111100011100010101001000010001111101011110110001010100101010110000011100110001001010001111110001110011101001010000101100010101111101001011111101011011011010011111010011111001100000101000110001111101000111001100000111111001100010101111001110100000101100000010011010010011011011110001000111110011011110011101101101110100001110000011110000001000011100110100000101010111010101000000111111111001000001110001101010001101101100110010101001100001100011100011001010101101011000010101001011011100110101111111001101111111110101001111001110100100001010000010100000001011101011100111111100001011101001000100001000011110010000011111000001000110111001010000111110000001000111011011011110001000111001111000001101110100110001001010100110110001110100000000010110100000100100001100101011011000111001001111110000011110110000011001000000111011111010000110000101010001101010011111011100010011101111101010001010000000101100011000010101010101001100000000101110011100001111111011100010110000010101001111010001010000000111100100100111101111000111111111110000101101101111011000100000111100111100011011000101000010000011000010001011010011000101101010000010000110001100001001010110101010001110001110111111011111111000000011111101110111110000110100110101100010001111000100101000101111010101011101100000101011001011000001010010011000110001111111011010111011000100000100000010000101001101000110110001111110011001000010110010111001010100111000011000011100011010100101011111111010010101100101110100111100111110000000111101011100101100111011011011000101010101011111010010101010011110011011011101010101110001101101111000110110001101011001101001100100010110101011011100010001100010001010011100101000011011011001000000010011110110101101100111110111100111000010000001101100100111101101101111001100001000010110011010001010011100001111100111110011010010111011010101110000000100111110110011001010011111001100001110000111101110111100000110010001111111110011000000100011111010101011101111001010111000111000110100011011111000111000111111100101010001100011010000101001101111110111101011110100110110110000001001001001000101101100110011001111100110110111001100001101100100100010111101011101001001110001100100000100110001101101000000001001000010010111001110010100100010101000011001011111111010011111010101101101111010011010111000011001000100001100110001000100011001010001100001011101011101001100101110110011010000000110100010111000001010011000011101110111011110101100110110001110100001000101011111001000010001010000001000110100101100100001001000010010010001001011000110110011000010100010000101010010110110011000010010011010010001001001001011000111011101101000111111000001110100110110110000001100010111001100111101001011110110001101011110110101100111011111011001100010000011011010010000101010101010101001111111110101110111001110000001001100110100000100000000001111010101100000101111110011111011000111100101011010110101011001100000110111110001101101100000110111001011111011100110001111101101010111110011110011101011001110001110000101101011000101110010110011101000010010100101010000100001000001110101001010111001110011001000011001100001010100111010110101100001111100001110001001110111001110110101111001010110101010010001010111010100000001001001110000100011111101101111001011011110111011011010101011101110000101110";

  printf("CHAR[] : %s\n", data);
  aevol::BitSet_SIMD** bitset = new aevol::BitSet_SIMD*[1024];

  for (int i = 0; i < 1024; i++)
    bitset[i] = new aevol::BitSet_SIMD(data, 5022);


  int8_t* dist_all = (int8_t*) _mm_malloc(1024 * 4096 * sizeof(int8_t),32);

  __assume_aligned(dist_all, 32);

  for (unsigned j = 0; j < 1024; j++) {
    for (unsigned i = 0; i < 4096; i++) {
      __declspec(align(32)) int8_t dist[22];
      unsigned cpt = i*5000+j;
      __assume(cpt%4096==0);
      __assume_aligned(dist_all,32);

      dist[0] = (CHECK_BIT(aevol::BITSET_PROM_SEQ_LEAD,0) == ARRAY_CHECK_BIT(bitset[i]->data_,i,sizeof(char))) ? 0 : 1;
      dist[1] = (CHECK_BIT(aevol::BITSET_PROM_SEQ_LEAD,1) == ARRAY_CHECK_BIT(bitset[i]->data_,i+1,sizeof(char))) ? 0 : 1;
      dist[2] = (CHECK_BIT(aevol::BITSET_PROM_SEQ_LEAD,2) == ARRAY_CHECK_BIT(bitset[i]->data_,i+2,sizeof(char))) ? 0 : 1;
      dist[3] = (CHECK_BIT(aevol::BITSET_PROM_SEQ_LEAD,3) == ARRAY_CHECK_BIT(bitset[i]->data_,i+3,sizeof(char))) ? 0 : 1;
      dist[4] = (CHECK_BIT(aevol::BITSET_PROM_SEQ_LEAD,4) == ARRAY_CHECK_BIT(bitset[i]->data_,i+4,sizeof(char))) ? 0 : 1;
      dist[5] = (CHECK_BIT(aevol::BITSET_PROM_SEQ_LEAD,5) == ARRAY_CHECK_BIT(bitset[i]->data_,i+5,sizeof(char))) ? 0 : 1;
      dist[6] = (CHECK_BIT(aevol::BITSET_PROM_SEQ_LEAD,6) == ARRAY_CHECK_BIT(bitset[i]->data_,i+6,sizeof(char))) ? 0 : 1;
      dist[7] = (CHECK_BIT(aevol::BITSET_PROM_SEQ_LEAD,7) == ARRAY_CHECK_BIT(bitset[i]->data_,i+7,sizeof(char))) ? 0 : 1;
      dist[8] = (CHECK_BIT(aevol::BITSET_PROM_SEQ_LEAD,8) == ARRAY_CHECK_BIT(bitset[i]->data_,i+8,sizeof(char))) ? 0 : 1;
      dist[9] = (CHECK_BIT(aevol::BITSET_PROM_SEQ_LEAD,9) == ARRAY_CHECK_BIT(bitset[i]->data_,i+9,sizeof(char))) ? 0 : 1;
      dist[10] = (CHECK_BIT(aevol::BITSET_PROM_SEQ_LEAD,10) == ARRAY_CHECK_BIT(bitset[i]->data_,i+10,sizeof(char))) ? 0 : 1;
      dist[11] = (CHECK_BIT(aevol::BITSET_PROM_SEQ_LEAD,11) == ARRAY_CHECK_BIT(bitset[i]->data_,i+11,sizeof(char))) ? 0 : 1;
      dist[12] = (CHECK_BIT(aevol::BITSET_PROM_SEQ_LEAD,12) == ARRAY_CHECK_BIT(bitset[i]->data_,i+12,sizeof(char))) ? 0 : 1;
      dist[13] = (CHECK_BIT(aevol::BITSET_PROM_SEQ_LEAD,13) == ARRAY_CHECK_BIT(bitset[i]->data_,i+13,sizeof(char))) ? 0 : 1;
      dist[14] = (CHECK_BIT(aevol::BITSET_PROM_SEQ_LEAD,14) == ARRAY_CHECK_BIT(bitset[i]->data_,i+14,sizeof(char))) ? 0 : 1;
      dist[15] = (CHECK_BIT(aevol::BITSET_PROM_SEQ_LEAD,15) == ARRAY_CHECK_BIT(bitset[i]->data_,i+15,sizeof(char))) ? 0 : 1;
      dist[16] = (CHECK_BIT(aevol::BITSET_PROM_SEQ_LEAD,16) == ARRAY_CHECK_BIT(bitset[i]->data_,i+16,sizeof(char))) ? 0 : 1;
      dist[17] = (CHECK_BIT(aevol::BITSET_PROM_SEQ_LEAD,17) == ARRAY_CHECK_BIT(bitset[i]->data_,i+17,sizeof(char))) ? 0 : 1;
      dist[18] = (CHECK_BIT(aevol::BITSET_PROM_SEQ_LEAD,18) == ARRAY_CHECK_BIT(bitset[i]->data_,i+18,sizeof(char))) ? 0 : 1;
      dist[19] = (CHECK_BIT(aevol::BITSET_PROM_SEQ_LEAD,19) == ARRAY_CHECK_BIT(bitset[i]->data_,i+19,sizeof(char))) ? 0 : 1;
      dist[20] = (CHECK_BIT(aevol::BITSET_PROM_SEQ_LEAD,20) == ARRAY_CHECK_BIT(bitset[i]->data_,i+20,sizeof(char))) ? 0 : 1;
      dist[21] = (CHECK_BIT(aevol::BITSET_PROM_SEQ_LEAD,21) == ARRAY_CHECK_BIT(bitset[i]->data_,i+21,sizeof(char))) ? 0 : 1;
      dist_all[i*5000+j] = dist[0] + dist[1]+ dist[2]+ dist[3]+ dist[4]+ dist[5]+ dist[6]
             + dist[7] + dist[8]+ dist[9]+ dist[10]+ dist[11]+ dist[12]+ dist[13]
             + dist[14]+ dist[15]+ dist[16]+ dist[17]+ dist[18]+ dist[19]+ dist[20]
             + dist[21];
    }
  }
}*/
/*
  aevol::BitSet_SIMD* clone_bitset = new aevol::BitSet_SIMD(bitset);
  clone_bitset->print();
  bitset->remove(2,4);
  bitset->print();

  char BITSET_PROTEIN_END_LEAD  = 0b1110;

  int motif_id = 7;
  int position = 1;

  for (int i = 0; i < 8; i++)
    if (BITSET_PROTEIN_END_LEAD & (1 << (i)))
      printf("1");
    else
      printf("0");
  printf("\n");

  int8_t val = ((BITSET_PROTEIN_END_LEAD >> motif_id) & 1);

  int8_t dist = ((BITSET_PROTEIN_END_LEAD >> motif_id) & 1) == bitset->get(position) ? 1 : 0;

  printf("%d %d %d\n",dist,val,bitset->get(position));
*//*

//  bitset->flip(1);
//  bitset->print();

  //bitset->append(clone_bitset);
  //bitset->insert_at(clone_bitset,3);
  //bitset->print();
}

*/
