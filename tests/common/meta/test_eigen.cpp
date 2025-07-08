#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "galileo/common/meta/eigen.hpp"
#include <Eigen/Dense>

using namespace galileo;

// Test fixture for common matrix/vector setups
class EigenTestFixture
{
public:
    // Test matrices of various sizes
    static constexpr int ROWS = 10;
    static constexpr int COLS = 8;

    galileo::Matrix<double, ROWS, COLS> test_matrix;
    galileo::Matrix<double, ROWS, 1> test_col_vector;
    galileo::Matrix<double, 1, COLS> test_row_vector;
    galileo::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> test_dynamic_matrix;
    galileo::Matrix<double, Eigen::Dynamic, 1> test_dynamic_col_vector;
    galileo::Matrix<double, 1, Eigen::Dynamic> test_dynamic_row_vector;

    EigenTestFixture()
    {
        // Initialize test matrices with predictable values
        test_matrix = galileo::Matrix<double, ROWS, COLS>::Random();
        test_col_vector = galileo::Matrix<double, ROWS, 1>::Random();
        test_row_vector = galileo::Matrix<double, 1, COLS>::Random();

        test_dynamic_matrix = galileo::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>::Random(ROWS, COLS);
        test_dynamic_col_vector = galileo::Matrix<double, Eigen::Dynamic, 1>::Random(ROWS, 1);
        test_dynamic_row_vector = galileo::Matrix<double, 1, Eigen::Dynamic>::Random(1, COLS);
    }
};

TEST_CASE("Matrix Storage Order Correction", "[eigen][matrix]")
{
    SECTION("Row vector should be row-major")
    {
        using RowVectorType = galileo::Matrix<double, 1, 5>;
        static_assert(static_cast<int>(RowVectorType::Options) & static_cast<int>(Eigen::RowMajor));
        static_assert(IsEigenRowVector<RowVectorType>);

        RowVectorType row_vec;
        REQUIRE(row_vec.rows() == 1);
        REQUIRE(row_vec.cols() == 5);
    }

    SECTION("Column vector should be column-major")
    {
        using ColVectorType = galileo::Matrix<double, 5, 1>;
        static_assert(!(static_cast<int>(ColVectorType::Options) & static_cast<int>(Eigen::RowMajor)));
        static_assert(IsEigenColVector<ColVectorType>);

        ColVectorType col_vec;
        REQUIRE(col_vec.rows() == 5);
        REQUIRE(col_vec.cols() == 1);
    }

    SECTION("Regular matrix preserves default storage order")
    {
        using MatrixType = galileo::Matrix<double, 3, 4>;
        static_assert(!(static_cast<int>(MatrixType::Options) & static_cast<int>(Eigen::RowMajor)));
        static_assert(IsEigenMatrix<MatrixType>);

        MatrixType mat;
        REQUIRE(mat.rows() == 3);
        REQUIRE(mat.cols() == 4);
    }

    SECTION("Dynamic row vector should be row-major")
    {
        using DynamicRowVectorType = galileo::Matrix<double, 1, Eigen::Dynamic>;
        static_assert(static_cast<int>(DynamicRowVectorType::Options) & static_cast<int>(Eigen::RowMajor));
        static_assert(IsEigenRowVector<DynamicRowVectorType>);

        DynamicRowVectorType dyn_row_vec(1, 6);
        REQUIRE(dyn_row_vec.rows() == 1);
        REQUIRE(dyn_row_vec.cols() == 6);
    }

    SECTION("Dynamic column vector should be column-major")
    {
        using DynamicColVectorType = galileo::Matrix<double, Eigen::Dynamic, 1>;
        static_assert(!(static_cast<int>(DynamicColVectorType::Options) & static_cast<int>(Eigen::RowMajor)));
        static_assert(IsEigenColVector<DynamicColVectorType>);

        DynamicColVectorType dyn_col_vec(6, 1);
        REQUIRE(dyn_col_vec.rows() == 6);
        REQUIRE(dyn_col_vec.cols() == 1);
    }
}

TEST_CASE_METHOD(EigenTestFixture, "Vector Segment Dispatchers", "[eigen][vector][segment]")
{
    SECTION("segment with compile-time size")
    {
        constexpr int segment_size = 3;
        constexpr int start_idx = 2;

        auto seg1 = segment<segment_size>(test_col_vector, start_idx);
        auto seg2 = segment<segment_size>(test_row_vector, start_idx);

        REQUIRE(seg1.size() == segment_size);
        REQUIRE(seg2.size() == segment_size);

        for (int i = 0; i < segment_size; ++i)
        {
            REQUIRE(seg1[i] == test_col_vector[start_idx + i]);
            REQUIRE(seg2[i] == test_row_vector[start_idx + i]);
        }
    }

    SECTION("segment with runtime size")
    {
        const int segment_size = 4;
        const int start_idx = 1;

        auto seg1 = segment(test_col_vector, start_idx, segment_size);
        auto seg2 = segment(test_row_vector, start_idx, segment_size);

        REQUIRE(seg1.size() == segment_size);
        REQUIRE(seg2.size() == segment_size);

        for (int i = 0; i < segment_size; ++i)
        {
            REQUIRE(seg1[i] == test_col_vector[start_idx + i]);
            REQUIRE(seg2[i] == test_row_vector[start_idx + i]);
        }
    }

    SECTION("segment with DimensionTpl<> (dynamic)")
    {
        DimensionTpl<> segment_dim(5);
        const int start_idx = 0;

        auto seg1 = segment(test_col_vector, start_idx, segment_dim);
        auto seg2 = segment(test_row_vector, start_idx, segment_dim);

        REQUIRE(seg1.size() == segment_dim.value());
        REQUIRE(seg2.size() == segment_dim.value());

        for (int i = 0; i < segment_dim.value(); ++i)
        {
            REQUIRE(seg1[i] == test_col_vector[start_idx + i]);
            REQUIRE(seg2[i] == test_row_vector[start_idx + i]);
        }
    }

    SECTION("segment with DimensionTpl<N> (fixed)")
    {
        constexpr DimensionTpl<4> segment_dim;
        const int start_idx = 3;

        auto seg1 = segment(test_col_vector, start_idx, segment_dim);
        auto seg2 = segment(test_row_vector, start_idx, segment_dim);

        REQUIRE(seg1.size() == segment_dim.value());
        REQUIRE(seg2.size() == segment_dim.value());

        for (int i = 0; i < segment_dim.value(); ++i)
        {
            REQUIRE(seg1[i] == test_col_vector[start_idx + i]);
            REQUIRE(seg2[i] == test_row_vector[start_idx + i]);
        }
    }

    SECTION("segment with dynamic vectors")
    {
        constexpr int segment_size = 3;
        const int start_idx = 1;

        auto seg1 = segment<segment_size>(test_dynamic_col_vector, start_idx);
        auto seg2 = segment<segment_size>(test_dynamic_row_vector, start_idx);

        REQUIRE(seg1.size() == segment_size);
        REQUIRE(seg2.size() == segment_size);

        for (int i = 0; i < segment_size; ++i)
        {
            REQUIRE(seg1[i] == test_dynamic_col_vector[start_idx + i]);
            REQUIRE(seg2[i] == test_dynamic_row_vector[start_idx + i]);
        }
    }
}

TEST_CASE_METHOD(EigenTestFixture, "Vector Head Dispatchers", "[eigen][vector][head]")
{
    SECTION("head with compile-time size")
    {
        constexpr int head_size = 3;

        auto head1 = head<head_size>(test_col_vector);
        auto head2 = head<head_size>(test_row_vector);

        REQUIRE(head1.size() == head_size);
        REQUIRE(head2.size() == head_size);

        for (int i = 0; i < head_size; ++i)
        {
            REQUIRE(head1[i] == test_col_vector[i]);
            REQUIRE(head2[i] == test_row_vector[i]);
        }
    }

    SECTION("head with runtime size")
    {
        const int head_size = 4;

        auto head1 = head(test_col_vector, head_size);
        auto head2 = head(test_row_vector, head_size);

        REQUIRE(head1.size() == head_size);
        REQUIRE(head2.size() == head_size);

        for (int i = 0; i < head_size; ++i)
        {
            REQUIRE(head1[i] == test_col_vector[i]);
            REQUIRE(head2[i] == test_row_vector[i]);
        }
    }

    SECTION("head with DimensionTpl<> (dynamic)")
    {
        DimensionTpl<> head_dim(5);

        auto head1 = head(test_col_vector, head_dim);
        auto head2 = head(test_row_vector, head_dim);

        REQUIRE(head1.size() == head_dim.value());
        REQUIRE(head2.size() == head_dim.value());

        for (int i = 0; i < head_dim.value(); ++i)
        {
            REQUIRE(head1[i] == test_col_vector[i]);
            REQUIRE(head2[i] == test_row_vector[i]);
        }
    }

    SECTION("head with DimensionTpl<N> (fixed)")
    {
        constexpr DimensionTpl<3> head_dim;

        auto head1 = head(test_col_vector, head_dim);
        auto head2 = head(test_row_vector, head_dim);

        REQUIRE(head1.size() == head_dim.value());
        REQUIRE(head2.size() == head_dim.value());

        for (int i = 0; i < head_dim.value(); ++i)
        {
            REQUIRE(head1[i] == test_col_vector[i]);
            REQUIRE(head2[i] == test_row_vector[i]);
        }
    }
}

TEST_CASE_METHOD(EigenTestFixture, "Vector Tail Dispatchers", "[eigen][vector][tail]")
{
    SECTION("tail with compile-time size")
    {
        constexpr int tail_size = 3;

        auto tail1 = tail<tail_size>(test_col_vector);
        auto tail2 = tail<tail_size>(test_row_vector);

        REQUIRE(tail1.size() == tail_size);
        REQUIRE(tail2.size() == tail_size);

        for (int i = 0; i < tail_size; ++i)
        {
            REQUIRE(tail1[i] == test_col_vector[ROWS - tail_size + i]);
            REQUIRE(tail2[i] == test_row_vector[COLS - tail_size + i]);
        }
    }

    SECTION("tail with runtime size")
    {
        const int tail_size = 4;

        auto tail1 = tail(test_col_vector, tail_size);
        auto tail2 = tail(test_row_vector, tail_size);

        REQUIRE(tail1.size() == tail_size);
        REQUIRE(tail2.size() == tail_size);

        for (int i = 0; i < tail_size; ++i)
        {
            REQUIRE(tail1[i] == test_col_vector[ROWS - tail_size + i]);
            REQUIRE(tail2[i] == test_row_vector[COLS - tail_size + i]);
        }
    }

    SECTION("tail with DimensionTpl<> (dynamic)")
    {
        DimensionTpl<> tail_dim(2);

        auto tail1 = tail(test_col_vector, tail_dim);
        auto tail2 = tail(test_row_vector, tail_dim);

        REQUIRE(tail1.size() == tail_dim.value());
        REQUIRE(tail2.size() == tail_dim.value());

        for (int i = 0; i < tail_dim.value(); ++i)
        {
            REQUIRE(tail1[i] == test_col_vector[ROWS - tail_dim.value() + i]);
            REQUIRE(tail2[i] == test_row_vector[COLS - tail_dim.value() + i]);
        }
    }

    SECTION("tail with DimensionTpl<N> (fixed)")
    {
        constexpr DimensionTpl<5> tail_dim;

        auto tail1 = tail(test_col_vector, tail_dim);
        auto tail2 = tail(test_row_vector, tail_dim);

        REQUIRE(tail1.size() == tail_dim.value());
        REQUIRE(tail2.size() == tail_dim.value());

        for (int i = 0; i < tail_dim.value(); ++i)
        {
            REQUIRE(tail1[i] == test_col_vector[ROWS - tail_dim.value() + i]);
            REQUIRE(tail2[i] == test_row_vector[COLS - tail_dim.value() + i]);
        }
    }
}

TEST_CASE_METHOD(EigenTestFixture, "Matrix Block Dispatchers", "[eigen][matrix][block]")
{
    SECTION("block with compile-time sizes")
    {
        constexpr int block_rows = 3;
        constexpr int block_cols = 4;
        const int start_row = 1;
        const int start_col = 2;

        auto block_result = block<block_rows, block_cols>(test_matrix, start_row, start_col);

        REQUIRE(block_result.rows() == block_rows);
        REQUIRE(block_result.cols() == block_cols);

        for (int i = 0; i < block_rows; ++i)
        {
            for (int j = 0; j < block_cols; ++j)
            {
                REQUIRE(block_result(i, j) == test_matrix(start_row + i, start_col + j));
            }
        }
    }

    SECTION("block with runtime sizes")
    {
        const int block_rows = 2;
        const int block_cols = 3;
        const int start_row = 0;
        const int start_col = 1;

        auto block_result = block(test_matrix, start_row, start_col, block_rows, block_cols);

        REQUIRE(block_result.rows() == block_rows);
        REQUIRE(block_result.cols() == block_cols);

        for (int i = 0; i < block_rows; ++i)
        {
            for (int j = 0; j < block_cols; ++j)
            {
                REQUIRE(block_result(i, j) == test_matrix(start_row + i, start_col + j));
            }
        }
    }

    SECTION("block with mixed compile-time and runtime sizes")
    {
        constexpr int block_rows = 4;
        const int block_cols = 2;
        const int start_row = 2;
        const int start_col = 3;

        auto block_result = block(test_matrix, start_row, start_col, block_rows, block_cols);

        REQUIRE(block_result.rows() == block_rows);
        REQUIRE(block_result.cols() == block_cols);

        for (int i = 0; i < block_rows; ++i)
        {
            for (int j = 0; j < block_cols; ++j)
            {
                REQUIRE(block_result(i, j) == test_matrix(start_row + i, start_col + j));
            }
        }
    }

    SECTION("block with DimensionTpl objects")
    {
        constexpr DimensionTpl<3> block_rows_dim;
        DimensionTpl<> block_cols_dim(2);
        const int start_row = 1;
        const int start_col = 0;

        auto block_result = block(test_matrix, start_row, start_col, block_rows_dim, block_cols_dim);

        REQUIRE(block_result.rows() == block_rows_dim.value());
        REQUIRE(block_result.cols() == block_cols_dim.value());

        for (int i = 0; i < block_rows_dim.value(); ++i)
        {
            for (int j = 0; j < block_cols_dim.value(); ++j)
            {
                REQUIRE(block_result(i, j) == test_matrix(start_row + i, start_col + j));
            }
        }
    }

    SECTION("block with dynamic matrix")
    {
        constexpr int block_rows = 2;
        constexpr int block_cols = 3;
        const int start_row = 3;
        const int start_col = 1;

        auto block_result = block<block_rows, block_cols>(test_dynamic_matrix, start_row, start_col);

        REQUIRE(block_result.rows() == block_rows);
        REQUIRE(block_result.cols() == block_cols);

        for (int i = 0; i < block_rows; ++i)
        {
            for (int j = 0; j < block_cols; ++j)
            {
                REQUIRE(block_result(i, j) == test_dynamic_matrix(start_row + i, start_col + j));
            }
        }
    }
}

TEST_CASE_METHOD(EigenTestFixture, "Matrix Corner Dispatchers", "[eigen][matrix][corner]")
{
    SECTION("topLeftCorner with compile-time sizes")
    {
        constexpr int corner_rows = 2;
        constexpr int corner_cols = 3;

        auto corner_result = topLeftCorner<corner_rows, corner_cols>(test_matrix);

        REQUIRE(corner_result.rows() == corner_rows);
        REQUIRE(corner_result.cols() == corner_cols);

        for (int i = 0; i < corner_rows; ++i)
        {
            for (int j = 0; j < corner_cols; ++j)
            {
                REQUIRE(corner_result(i, j) == test_matrix(i, j));
            }
        }
    }

    SECTION("topLeftCorner with runtime sizes")
    {
        const int corner_rows = 3;
        const int corner_cols = 2;

        auto corner_result = topLeftCorner(test_matrix, corner_rows, corner_cols);

        REQUIRE(corner_result.rows() == corner_rows);
        REQUIRE(corner_result.cols() == corner_cols);

        for (int i = 0; i < corner_rows; ++i)
        {
            for (int j = 0; j < corner_cols; ++j)
            {
                REQUIRE(corner_result(i, j) == test_matrix(i, j));
            }
        }
    }

    SECTION("topRightCorner with DimensionTpl objects")
    {
        constexpr DimensionTpl<2> corner_rows_dim;
        DimensionTpl<> corner_cols_dim(4);

        auto corner_result = topRightCorner(test_matrix, corner_rows_dim, corner_cols_dim);

        REQUIRE(corner_result.rows() == corner_rows_dim.value());
        REQUIRE(corner_result.cols() == corner_cols_dim.value());

        for (int i = 0; i < corner_rows_dim.value(); ++i)
        {
            for (int j = 0; j < corner_cols_dim.value(); ++j)
            {
                REQUIRE(corner_result(i, j) == test_matrix(i, COLS - corner_cols_dim.value() + j));
            }
        }
    }

    SECTION("bottomLeftCorner with mixed sizes")
    {
        constexpr int corner_rows = 3;
        const int corner_cols = 2;

        auto corner_result = bottomLeftCorner(test_matrix, corner_rows, corner_cols);

        REQUIRE(corner_result.rows() == corner_rows);
        REQUIRE(corner_result.cols() == corner_cols);

        for (int i = 0; i < corner_rows; ++i)
        {
            for (int j = 0; j < corner_cols; ++j)
            {
                REQUIRE(corner_result(i, j) == test_matrix(ROWS - corner_rows + i, j));
            }
        }
    }

    SECTION("bottomRightCorner with compile-time sizes")
    {
        constexpr int corner_rows = 4;
        constexpr int corner_cols = 3;

        auto corner_result = bottomRightCorner<corner_rows, corner_cols>(test_matrix);

        REQUIRE(corner_result.rows() == corner_rows);
        REQUIRE(corner_result.cols() == corner_cols);

        for (int i = 0; i < corner_rows; ++i)
        {
            for (int j = 0; j < corner_cols; ++j)
            {
                REQUIRE(corner_result(i, j) == test_matrix(ROWS - corner_rows + i, COLS - corner_cols + j));
            }
        }
    }
}

TEST_CASE_METHOD(EigenTestFixture, "Matrix Row/Column Dispatchers", "[eigen][matrix][rows_cols]")
{
    SECTION("topRows with compile-time size")
    {
        constexpr int num_rows = 3;

        auto rows_result = topRows<num_rows>(test_matrix);

        REQUIRE(rows_result.rows() == num_rows);
        REQUIRE(rows_result.cols() == COLS);

        for (int i = 0; i < num_rows; ++i)
        {
            for (int j = 0; j < COLS; ++j)
            {
                REQUIRE(rows_result(i, j) == test_matrix(i, j));
            }
        }
    }

    SECTION("topRows with runtime size")
    {
        const int num_rows = 4;

        auto rows_result = topRows(test_matrix, num_rows);

        REQUIRE(rows_result.rows() == num_rows);
        REQUIRE(rows_result.cols() == COLS);

        for (int i = 0; i < num_rows; ++i)
        {
            for (int j = 0; j < COLS; ++j)
            {
                REQUIRE(rows_result(i, j) == test_matrix(i, j));
            }
        }
    }

    SECTION("bottomRows with DimensionTpl object")
    {
        constexpr DimensionTpl<2> num_rows_dim;

        auto rows_result = bottomRows(test_matrix, num_rows_dim);

        REQUIRE(rows_result.rows() == num_rows_dim.value());
        REQUIRE(rows_result.cols() == COLS);

        for (int i = 0; i < num_rows_dim.value(); ++i)
        {
            for (int j = 0; j < COLS; ++j)
            {
                REQUIRE(rows_result(i, j) == test_matrix(ROWS - num_rows_dim.value() + i, j));
            }
        }
    }

    SECTION("leftCols with compile-time size")
    {
        constexpr int num_cols = 3;

        auto cols_result = leftCols<num_cols>(test_matrix);

        REQUIRE(cols_result.rows() == ROWS);
        REQUIRE(cols_result.cols() == num_cols);

        for (int i = 0; i < ROWS; ++i)
        {
            for (int j = 0; j < num_cols; ++j)
            {
                REQUIRE(cols_result(i, j) == test_matrix(i, j));
            }
        }
    }

    SECTION("rightCols with dynamic DimensionTpl")
    {
        DimensionTpl<> num_cols_dim(4);

        auto cols_result = rightCols(test_matrix, num_cols_dim);

        REQUIRE(cols_result.rows() == ROWS);
        REQUIRE(cols_result.cols() == num_cols_dim.value());

        for (int i = 0; i < ROWS; ++i)
        {
            for (int j = 0; j < num_cols_dim.value(); ++j)
            {
                REQUIRE(cols_result(i, j) == test_matrix(i, COLS - num_cols_dim.value() + j));
            }
        }
    }

    SECTION("Mixed operations with dynamic matrix")
    {
        const int num_rows = 2;
        const int num_cols = 3;

        auto top_rows = topRows(test_dynamic_matrix, num_rows);
        auto left_cols = leftCols(test_dynamic_matrix, num_cols);

        REQUIRE(top_rows.rows() == num_rows);
        REQUIRE(top_rows.cols() == COLS);
        REQUIRE(left_cols.rows() == ROWS);
        REQUIRE(left_cols.cols() == num_cols);

        for (int i = 0; i < num_rows; ++i)
        {
            for (int j = 0; j < COLS; ++j)
            {
                REQUIRE(top_rows(i, j) == test_dynamic_matrix(i, j));
            }
        }

        for (int i = 0; i < ROWS; ++i)
        {
            for (int j = 0; j < num_cols; ++j)
            {
                REQUIRE(left_cols(i, j) == test_dynamic_matrix(i, j));
            }
        }
    }
}

TEST_CASE("Edge Cases and Comprehensive Dispatch Testing", "[eigen][edge_cases]")
{
    SECTION("Zero-sized operations")
    {
        galileo::Matrix<double, 5, 5> test_mat = galileo::Matrix<double, 5, 5>::Random();

        // Zero-size head/tail
        auto zero_head = head<0>(test_mat.col(0));
        auto zero_tail = tail<0>(test_mat.col(0));

        REQUIRE(zero_head.size() == 0);
        REQUIRE(zero_tail.size() == 0);

        // Zero-size block
        auto zero_block = block<0, 0>(test_mat, 0, 0);
        REQUIRE(zero_block.rows() == 0);
        REQUIRE(zero_block.cols() == 0);
    }

    SECTION("Single element operations")
    {
        galileo::Matrix<double, 5, 5> test_mat = galileo::Matrix<double, 5, 5>::Random();

        auto single_head = head<1>(test_mat.col(0));
        auto single_tail = tail<1>(test_mat.col(0));
        auto single_segment = segment<1>(test_mat.col(0), 2);

        REQUIRE(single_head.size() == 1);
        REQUIRE(single_tail.size() == 1);
        REQUIRE(single_segment.size() == 1);

        REQUIRE(single_head[0] == test_mat(0, 0));
        REQUIRE(single_tail[0] == test_mat(4, 0));
        REQUIRE(single_segment[0] == test_mat(2, 0));
    }

    SECTION("Full-size operations")
    {
        constexpr int SIZE = 6;
        galileo::Matrix<double, SIZE, SIZE> test_mat = galileo::Matrix<double, SIZE, SIZE>::Random();

        auto full_head = head<SIZE>(test_mat.col(0));
        auto full_tail = tail<SIZE>(test_mat.col(0));
        auto full_segment = segment<SIZE>(test_mat.col(0), 0);

        REQUIRE(full_head.size() == SIZE);
        REQUIRE(full_tail.size() == SIZE);
        REQUIRE(full_segment.size() == SIZE);

        for (int i = 0; i < SIZE; ++i)
        {
            REQUIRE(full_head[i] == test_mat(i, 0));
            REQUIRE(full_tail[i] == test_mat(i, 0));
            REQUIRE(full_segment[i] == test_mat(i, 0));
        }
    }

    SECTION("Complex chained operations")
    {
        galileo::Matrix<double, 8, 8> test_mat = galileo::Matrix<double, 8, 8>::Random();

        // Chain multiple operations
        auto complex_result = topRows<4>(leftCols<6>(test_mat));
        auto block_result = block<4, 6>(test_mat, 0, 0);

        REQUIRE(complex_result.rows() == 4);
        REQUIRE(complex_result.cols() == 6);
        REQUIRE(block_result.rows() == 4);
        REQUIRE(block_result.cols() == 6);

        for (int i = 0; i < 4; ++i)
        {
            for (int j = 0; j < 6; ++j)
            {
                REQUIRE(complex_result(i, j) == test_mat(i, j));
                REQUIRE(block_result(i, j) == test_mat(i, j));
            }
        }
    }

    SECTION("All dispatcher combinations with different types")
    {
        constexpr int TEST_SIZE = 6;
        galileo::Matrix<double, TEST_SIZE, TEST_SIZE> test_mat = galileo::Matrix<double, TEST_SIZE, TEST_SIZE>::Random();

        // Test all combinations of compile-time constants, runtime ints, and DimensionTpl objects
        constexpr int compile_time_size = 2;
        const int runtime_size = 3;
        constexpr DimensionTpl<2> fixed_dim;
        DimensionTpl<> dynamic_dim(3);

        // Test segment with all combinations
        auto seg1 = segment<compile_time_size>(test_mat.col(0), 0);
        auto seg2 = segment(test_mat.col(0), 1, runtime_size);
        auto seg3 = segment(test_mat.col(0), 2, fixed_dim);
        auto seg4 = segment(test_mat.col(0), 3, dynamic_dim);

        REQUIRE(seg1.size() == compile_time_size);
        REQUIRE(seg2.size() == runtime_size);
        REQUIRE(seg3.size() == fixed_dim.value());
        REQUIRE(seg4.size() == dynamic_dim.value());

        // Test block with all combinations
        auto block1 = block<compile_time_size, compile_time_size>(test_mat, 0, 0);
        auto block2 = block(test_mat, 0, 1, runtime_size, runtime_size);
        auto block3 = block(test_mat, 1, 0, fixed_dim, fixed_dim);
        auto block4 = block(test_mat, 1, 1, dynamic_dim, dynamic_dim);
        auto block5 = block(test_mat, 2, 0, compile_time_size, runtime_size);
        auto block6 = block(test_mat, 0, 2, fixed_dim, dynamic_dim);

        REQUIRE(block1.rows() == compile_time_size);
        REQUIRE(block1.cols() == compile_time_size);
        REQUIRE(block2.rows() == runtime_size);
        REQUIRE(block2.cols() == runtime_size);
        REQUIRE(block3.rows() == fixed_dim.value());
        REQUIRE(block3.cols() == fixed_dim.value());
        REQUIRE(block4.rows() == dynamic_dim.value());
        REQUIRE(block4.cols() == dynamic_dim.value());
        REQUIRE(block5.rows() == compile_time_size);
        REQUIRE(block5.cols() == runtime_size);
        REQUIRE(block6.rows() == fixed_dim.value());
        REQUIRE(block6.cols() == dynamic_dim.value());
    }
}

TEST_CASE("Type Safety and Concept Validation", "[eigen][concepts]")
{
    SECTION("Concept validation")
    {
        using ColVectorType = galileo::Matrix<double, 5, 1>;
        using RowVectorType = galileo::Matrix<double, 1, 5>;
        using MatrixType = galileo::Matrix<double, 5, 5>;

        static_assert(IsEigenColVector<ColVectorType>);
        static_assert(!IsEigenRowVector<ColVectorType>);
        static_assert(IsEigenVector<ColVectorType>);
        static_assert(!IsEigenMatrix<ColVectorType>);

        static_assert(IsEigenRowVector<RowVectorType>);
        static_assert(!IsEigenColVector<RowVectorType>);
        static_assert(IsEigenVector<RowVectorType>);
        static_assert(!IsEigenMatrix<RowVectorType>);

        static_assert(!IsEigenRowVector<MatrixType>);
        static_assert(!IsEigenColVector<MatrixType>);
        static_assert(!IsEigenVector<MatrixType>);
        static_assert(IsEigenMatrix<MatrixType>);
    }

    SECTION("Extract compile time value validation")
    {
        static_assert(extract_compile_time_value<42>::Value == 42);
        static_assert(extract_compile_time_value<15>::Value == 15);
        static_assert(extract_compile_time_value<0>::Value == 0);

        // Note: extract_compile_time_value with DimensionTpl objects as template parameters
        // requires constexpr DimensionTpl values, but DimensionTpl objects cannot be used
        // as template non-type parameters in this context. This is by design.
        constexpr DimensionTpl<15> test_dim_15;
        constexpr DimensionTpl<0> test_dim_0;
        REQUIRE(test_dim_15.Value == 15);
        REQUIRE(test_dim_0.Value == 0);
    }

    SECTION("DimensionTpl arithmetic type correctness")
    {
        constexpr DimensionTpl<5> a;
        constexpr DimensionTpl<3> b;
        DimensionTpl<> c(7);

        static_assert(decltype(a + b)::Value == 8);
        static_assert(decltype(a + b)::IsFixed);
        static_assert(decltype(a + c)::IsDynamic);
        static_assert(decltype(c + b)::IsDynamic);
        static_assert(decltype(c + c)::IsDynamic);
    }
}
