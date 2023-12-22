classdef BlockMatrix < matlab.mixin.indexing.RedefinesParen
    %BLOCKMATRIX This class represents a block matrix.
    %   A block matrix is a matrix that is partitioned into smaller matrices,
    %   or blocks, which can be processed individually.
    
    properties
        M; % Number of rows in each block.
        N; % Number of columns in each block.
        A; % The actual block matrix stored as a 2D array.
    end

    methods
        function self = BlockMatrix(DimRow, DimCol)
            %BLOCKMATRIX Constructor for BlockMatrix class.
            %   Creates a new BlockMatrix object with specified dimensions for each block.
            %   DimRow: Number of rows in each block.
            %   DimCol: Number of columns in each block.
            self.M = DimRow;
            self.N = DimCol;
        end
    end
    
    methods (Access = protected)
        
        function varargout = parenReference(obj, indexOp)
            %PARENREFERENCE Overloaded method to reference block matrix elements.
            %   Retrieves elements from the block matrix based on expanded indices.
            %   indexOp: Cell array of indices to reference.

            indexOp = obj.expand(indexOp);
            varargout{1} = obj.A(indexOp{1},indexOp{2});
        end

        function obj = parenAssign(obj,indexOp,varargin)
            %PARENASSIGN Overloaded method to assign values to block matrix elements.
            %   Assigns values to the block matrix based on expanded indices.
            %   indexOp: Cell array of indices for assignment.
            %   varargin: Values to assign to the specified indices.

            indexOp = obj.expand(indexOp);
            % Ensure object instance is the first argument of call.
            if isempty(obj)
                obj.A = varargin{1};
                return;
            end
            assert(nargin==3);
            rhs = varargin{1};
            obj.A(indexOp{1},indexOp{2}) = rhs;
        end

        function n = parenListLength(obj,indexOp,ctx)
            %PARENLISTLENGTH Overloaded method to determine the number of elements to return.
            %   Determines the number of elements to return based on expanded indices.
            %   indexOp: Cell array of indices to determine the list length.
            %   ctx: Context of the indexing operation.
            indexOp = obj.expand(indexOp);
            if numel(indexOp) <= 2
                n = 1;
                return;
            end
            containedObj = obj.(indexOp(1:2));
            n = listLength(containedObj,indexOp(3:end),ctx);
        end

        function obj = parenDelete(obj,indexOp)
            %PARENDELETE Overloaded method to delete block matrix elements.
            %   Deletes elements from the block matrix based on expanded indices.
            %   indexOp: Cell array of indices to delete.
            indexOp = obj.expand(indexOp);
            obj.A(indexOp{1},indexOp{2}) = [];
        end
    end

    methods (Access=private)
        function indexOpMod = expand(self, indexOp)
            %EXPAND Helper method to expand block indices to matrix indices.
            %   Expands block indices to the corresponding matrix indices for referencing or assignment.
            %   indexOp: Original block indices to expand.

            idxRow = [];
            idxCol = [];
            for e = indexOp.Indices{1}
                el = self.M*(e-1)+1:self.M*e;
                idxRow = [idxRow,el];
            end
            for e = indexOp.Indices{2}
                el = self.N*(e-1)+1:self.N*e;
                idxCol = [idxCol,el];
            end
            indexOpMod = {idxRow, idxCol};
        end
    end

    methods (Access=public)
        function out = value(obj)
            %VALUE Method to get the value of the block matrix.
            %   Returns the actual block matrix stored as a 2D array.
            out = obj.A;
        end
        
        function out = sum(obj)
            %SUM Method to calculate the sum of all elements in the block matrix.
            %   Returns the sum of all elements in the block matrix.
            out = sum(obj.A,"all");
        end
        
        function out = cat(dim,varargin)
            %CAT Overloaded method to concatenate block matrices.
            %   Concatenates block matrices along the specified dimension.
            %   dim: Dimension along which to concatenate.
            %   varargin: Block matrices to concatenate.

            numCatArrays = nargin-1;
            newArgs = cell(numCatArrays,1);
            for ix = 1:numCatArrays
                if isa(varargin{ix},'utils.BlockMatrix')
                    newArgs{ix} = varargin{ix}.A;
                else
                    newArgs{ix} = varargin{ix};
                end
            end
            out = utils.BlockMatrix(varargin{1}.M,varargin{1}.N);
            out.A = cat(dim,newArgs{:});
        end

        function varargout = size(obj,varargin)
            %SIZE Overloaded method to get the size of the block matrix.
            %   Returns the size of the block matrix along specified dimensions.
            %   varargin: Dimensions to get the size of.
            [varargout{1:nargout}] = size(obj.A,varargin{:});
        end
    end

    methods (Static, Access=public)
        function obj = empty()
            %EMPTY Static method to create an empty BlockMatrix object.
            %   Returns an empty BlockMatrix object with zero dimensions for blocks.
            obj = utils.BlockMatrix(0,0);
        end
    end
end

