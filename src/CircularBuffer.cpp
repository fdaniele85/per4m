//
// Created by daniele on 14/02/25.
//

#include <per4m/CircularBuffer.h>

namespace ffp {
    CircularBuffer::CircularBuffer(const size_type capacity) : capacity_(capacity), data_(new double[capacity]) {}

    CircularBuffer::~CircularBuffer() { delete[] data_; }

    void CircularBuffer::add(const double value) {
        last_ = last_ % capacity_;
        data_[last_] = value;
        ++last_;
    }

    double *CircularBuffer::begin() { return data_; }
    double *CircularBuffer::end() { return data_ + capacity_; }
    const double *CircularBuffer::cbegin() const { return data_; }
    const double *CircularBuffer::cend() const { return data_ + capacity_; }
    const double *CircularBuffer::begin() const { return cbegin(); }
    const double *CircularBuffer::end() const { return cend(); }

    bool CircularBuffer::full() const { return last_ == capacity_; }

    double CircularBuffer::operator[](const size_type index) const { return data_[index]; }
} // namespace ffp