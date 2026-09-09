function varargout = trunc_mat_min(varargin)
%Truncate multidimensional matrices of same number dimensions to minimum size

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Sample usage
%[A, B] = trunc_mat_min(randn([3 4 5]), randn([2 8 3] ))

min_size = size(varargin{1});
for n = 2:nargin
    min_size = min(min_size, size(varargin{n}));
end

N_dims = numel(min_size);
idx_cell = cell(1, N_dims);
for i = 1:N_dims
    idx_cell{i} = 1:min_size(i);
end

varargout = cell(1, nargin);
for n = 1:nargin
    varargout{n} = varargin{n}(idx_cell{:});
end