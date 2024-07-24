% order the rows of a matrix A based on the values in column n
function  A_ordered = order_matrix_rows_from_column(A, column_index)
% A is the input matrix whose rows needed to be re-ordered
% column_index is the column index that will be refered to re-order
    [~, index] = sort(A(:, column_index));
    A_ordered = A(index, :);
end