function parts = find_tree_partitions_K6()
% FIND_TREE_PARTITIONS_K6
% Finds all ordered tree partitions of K6
% Automatically saves results to tree_partitions_K6.txt

    n = 6;
    [edgeList, edgeKeyToIndex] = make_edge_indexing(n);

    fprintf('Generating all spanning trees of K6...\n');
    [treeEdges, treeMask] = all_spanning_trees_Kn(n, edgeKeyToIndex);
    numTrees = size(treeEdges,1);
    fprintf('Total spanning trees: %d\n', numTrees);

    fullMask = uint32(2^15 - 1);
    parts = zeros(0, 15);

    fprintf('Searching for tree partitions...\n');

    for a = 1:numTrees
        maskA = treeMask(a);
        remAfterA = bitand(fullMask, bitcmp(maskA,'uint32'));

        for b = 1:numTrees
            maskB = treeMask(b);

            if bitand(maskA, maskB) ~= 0
                continue;
            end

            remAfterB = bitand(remAfterA, bitcmp(maskB,'uint32'));

            % Check if remaining edges form a spanning tree
            for c = 1:numTrees
                if treeMask(c) == remAfterB
                    A = sort(treeEdges(a,:));
                    B = sort(treeEdges(b,:));
                    C = sort(treeEdges(c,:));

                    row = [A B C];
                    parts(end+1,:) = row; %#ok<AGROW>
                end
            end
        end

        if mod(a,100)==0
            fprintf('Processed %d / %d trees...\n', a, numTrees);
        end
    end

    fprintf('\nTotal tree partitions found: %d\n', size(parts,1));

    % SAVE FILE
    writematrix(parts,'tree_partitions_K6.txt','Delimiter','space');
    fprintf('Saved to tree_partitions_K6.txt\n');

end

% ===============================
% Helper Functions (MUST BE BELOW)
% ===============================

function [edgeList, edgeKeyToIndex] = make_edge_indexing(n)

    edgeList = zeros(n*(n-1)/2, 2);
    idx = 0;
    for u = 1:n
        for v = u+1:n
            idx = idx + 1;
            edgeList(idx,:) = [u v];
        end
    end

    edgeKeyToIndex = containers.Map('KeyType','char','ValueType','int32');
    for i = 1:size(edgeList,1)
        u = edgeList(i,1); 
        v = edgeList(i,2);
        edgeKeyToIndex(sprintf('%d-%d',u,v)) = int32(i);
    end
end

function [treeEdges, treeMask] = all_spanning_trees_Kn(n, edgeKeyToIndex)

    numTrees = n^(n-2);  % Cayley formula
    treeEdges = zeros(numTrees, n-1);
    treeMask  = zeros(numTrees, 1, 'uint32');

    L = n-2;
    seq = zeros(1,L);

    for t = 1:numTrees

        % Convert t-1 to base n
        x = t-1;
        for i = 1:L
            seq(i) = mod(x,n) + 1;
            x = floor(x/n);
        end

        edgesUV = prufer_to_tree_edges(seq, n);

        eIdx = zeros(1,n-1);
        mask = uint32(0);

        for k = 1:n-1
            u = edgesUV(k,1);
            v = edgesUV(k,2);
            if u>v
                temp=u; u=v; v=temp;
            end
            key = sprintf('%d-%d',u,v);
            ei = edgeKeyToIndex(key);

            eIdx(k) = ei;
            mask = bitor(mask, bitshift(uint32(1), ei-1));
        end

        treeEdges(t,:) = sort(eIdx);
        treeMask(t) = mask;
    end
end

function edges = prufer_to_tree_edges(seq, n)

    L = n-2;
    deg = ones(1,n);

    for i = 1:L
        deg(seq(i)) = deg(seq(i)) + 1;
    end

    edges = zeros(n-1,2);

    for i = 1:L
        leaf = find(deg==1,1,'first');
        edges(i,:) = [leaf, seq(i)];
        deg(leaf) = deg(leaf)-1;
        deg(seq(i)) = deg(seq(i))-1;
    end

    last = find(deg==1);
    edges(n-1,:) = [last(1), last(2)];
end