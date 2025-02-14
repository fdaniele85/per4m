//
// Created by daniele on 14/02/25.
//

#include "CircularBuffer.h"

namespace ffp {
    CircularBuffer::CircularBuffer(int capacity) : capacity_(capacity), data_(new double[capacity]) {}

    CircularBuffer::~CircularBuffer() { delete[] data_; }

    void CircularBuffer::add(double value) {
        last_ = last_ % capacity_;
        data_[last_] = value;
        ++last_;
    }

    double *CircularBuffer::begin() { return data_; }
    double *CircularBuffer::end() { return data_ + capacity_; }
    const double *CircularBuffer::cbegin() const { return data_; }
    const double * CircularBuffer::cend() const { return data_ + capacity_; }

    bool CircularBuffer::full() const { return last_ == capacity_; }

    double CircularBuffer::operator[](const unsigned int index) const { return data_[index]; }
} // namespace ffp