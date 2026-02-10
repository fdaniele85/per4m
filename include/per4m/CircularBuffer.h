//
// Created by daniele on 14/02/25.
//

#pragma once

#include <cstddef>

namespace ffp {

    /**
     * \brief A circular buffer class for storing double values.
     */
    class CircularBuffer {
    public:
        using size_type = std::size_t;
        /**
         * \brief Constructs a CircularBuffer with the specified capacity.
         *
         * \param capacity The maximum number of elements the buffer can hold.
         */
        explicit CircularBuffer(size_type capacity);

        /**
         * \brief Destroys the CircularBuffer object.
         */
        ~CircularBuffer();

        /**
         * \brief Adds a value to the buffer.
         *
         * \param value The value to be added to the buffer.
         */
        void add(double value);

        /**
         * \brief Returns a pointer to the beginning of the buffer.
         *
         * \return A pointer to the first element in the buffer.
         */
        double *begin();

        /**
         * \brief Returns a pointer to the end of the buffer.
         *
         * \return A pointer to the element following the last element in the buffer.
         */
        double *end();

        /**
         * \brief Returns a constant pointer to the beginning of the buffer.
         *
         * \return A constant pointer to the first element in the buffer.
         */
        [[nodiscard]] const double *cbegin() const;
        [[nodiscard]] const double *begin() const;

        /**
         * \brief Returns a constant pointer to the end of the buffer.
         *
         * \return A constant pointer to the element following the last element in the buffer.
         */
        [[nodiscard]] const double *cend() const;
        [[nodiscard]] const double *end() const;

        /**
         * \brief Checks if the buffer is full.
         *
         * \return True if the buffer is full, false otherwise.
         */
        [[nodiscard]] bool full() const;

        /**
         * \brief Accesses the element at the specified index.
         *
         * \param index The index of the element to access.
         * \return The value of the element at the specified index.
         */
        [[nodiscard]] double operator[](size_type index) const;

    private:
        size_type capacity_{0}; ///< The maximum capacity of the buffer.
        size_type last_{0};     ///< The index of the last element in the buffer.
        double *data_;    ///< Pointer to the buffer's data array.
    };

} // namespace ffp
