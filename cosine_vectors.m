function cosine_val = cosine_vectors(vector_1, vector_2)
% This function is to compute a cosin between two vectors
if (size(vector_1,1) > 1)&&(size(vector_2,1) > 1)
    warning('current version no support');
    return;
end

cosine_val = sum(bsxfun(@times, vector_1, vector_2),2);