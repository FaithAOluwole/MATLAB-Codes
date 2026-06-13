function parts = find_tree_partitions_K6_ordered()
% FIND_TREE_PARTITIONS_K6_ORDERED
% Finds UNORDERED tree partitions of K6 under the rule:
%  - each block is sorted increasingly
%  - block first elements satisfy i1 < j1 < k1
% Automatically saves results to tree_partitions_K6_ordered.txt

    n = 6;
    [~, edgeKeyToIndex] = make_edge_indexing(n);

    fprintf('Generating all spanning trees of K6...\n');
    [treeEdges, treeMask] = all_spanning_trees_Kn(n, edgeKeyToIndex);
    numTrees = size(treeEdges,1);
    fprintf('Total spanning trees: %d\n', numTrees);

    fullMask = uint32(2^15 - 1);

    % Use a map to avoid duplicates (canonical row string -> true)
    seen = containers.Map('KeyType','char','ValueType','logical');
    parts = zeros(0, 15);

    fprintf('Searching for ORDERED tree partitions...\n');

    % Precompute a fast mask->index lookup (for finding c quickly)
    maskToIndex = containers.Map('KeyType','uint32','ValueType','int32');
    for t = 1:numTrees
        maskToIndex(treeMask(t)) = int32(t);
    end

    for a = 1:numTrees
        maskA = treeMask(a);
        remAfterA = bitand(fullMask, bitcmp(maskA,'uint32'));

        for b = 1:numTrees
            maskB = treeMask(b);

            if bitand(maskA, maskB) ~= 0
                continue;
            end

            remAfterB = bitand(remAfterA, bitcmp(maskB,'uint32'));

            % Find c such that treeMask(c) == remAfterB (if it exists)
            if ~isKey(maskToIndex, remAfterB)
                continue;
            end
            c = maskToIndex(remAfterB);

            % Build blocks and canonicalize ordering by first element
            A = sort(treeEdges(a,:));
            B = sort(treeEdges(b,:));
            C = sort(treeEdges(c,:));

            row = canonicalize_blocks_3(A, B, C); % enforces i1<j1<k1

            key = row_key(row);
            if ~isKey(seen, key)
                seen(key) = true;
                parts(end+1,:) = row; %#ok<AGROW>
            end
        end

        if mod(a,100)==0
            fprintf('Processed %d / %d trees...\n', a, numTrees);
        end
    end

    fprintf('\nTotal ORDERED tree partitions found: %d\n', size(parts,1));

    % SAVE FILE
    outname = 'tree_partitions_K6_ordered.txt';
    writematrix(parts, outname, 'Delimiter', 'space');
    fprintf('Saved to %s\n', outname);
end

% ===============================
% Helper Functions (MUST BE BELOW)
% ===============================

function row = canonicalize_blocks_3(A, B, C)
% Sort blocks by their first entry, so i1<j1<k1
    firsts = [A(1), B(1), C(1)];
    [~, idx] = sort(firsts);

    blocks = {A, B, C};
    row = [blocks{idx(1)}, blocks{idx(2)}, blocks{idx(3)}];
end

function k = row_key(row)
% Unique string key for deduping
    k = sprintf('%d,', row); % trailing comma is fine
end

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
            if u>v, temp=u; u=v; v=temp; end

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