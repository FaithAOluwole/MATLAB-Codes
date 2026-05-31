# MATLAB-Codes

In this section, we give the MATLAB code for $H_3$. 
The following code gives you the unordered cycle-free partitions of $K_6$

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


Next, we have the code that describes the group action as seen in Lemma \ref{lemma 2}

function w = triangle_move_tree_partition_K6_ordered(v, x, y, z)
%TRIANGLE_MOVE_TREE_PARTITION_K6_ORDERED
% Same triangle move as your unordered version, BUT outputs the result
% in the ORDERED (canonical) format:
%   - sort within each block
%   - then reorder blocks so that i1 < j1 < k1 (first entries strictly increasing)

    % ---- validate inputs ----
    if ~isvector(v) || numel(v) ~= 15
        error('v must be a 1x15 (or 15x1) vector.');
    end
    v = v(:).'; % row

    if ~(isscalar(x) && isscalar(y) && isscalar(z))
        error('x, y, z must be scalars.');
    end
    if ~(x < y && y < z)
        error('Need x<y<z.');
    end
    if any([x y z] < 1) || any([x y z] > 6)
        error('For K6, vertices must satisfy 1 <= x<y<z <= 6.');
    end

    if ~isequal(sort(v), 1:15)
        error('v must contain each label 1..15 exactly once.');
    end

    if ~is_tree_partition_K6(v)
        error('Input v is not a tree partition of K6.');
    end

    % ---- triangle edges ----
    n = 6;
    triEdges = [edge_label(x,y,n), edge_label(x,z,n), edge_label(y,z,n)];

    % ---- split v into 3 trees (blocks) ----
    T1 = v(1:5);
    T2 = v(6:10);
    T3 = v(11:15);

    % mem(t) in {1,2,3} tells which tree currently contains triEdges(t)
    mem = zeros(1,3);
    for t = 1:3
        e = triEdges(t);
        if ismember(e,T1)
            mem(t)=1;
        elseif ismember(e,T2)
            mem(t)=2;
        else
            mem(t)=3;
        end
    end

    permsIdx = perms(1:3);  % 6 permutations of positions
    found = false;
    w = [];

    edge_verts = edge_vertex_table(n);

    for r = 1:size(permsIdx,1)
        newMem = mem(permsIdx(r,:));

        % Must differ on at least TWO of the three triangle edges
        if sum(newMem ~= mem) < 2
            continue;
        end

        % Remove all triangle edges first
        U1 = T1; U2 = T2; U3 = T3;
        for t = 1:3
            e = triEdges(t);
            U1(U1==e) = [];
            U2(U2==e) = [];
            U3(U3==e) = [];
        end

        % Add them back according to newMem
        for t = 1:3
            e = triEdges(t);
            if newMem(t) == 1
                U1 = [U1, e];
            elseif newMem(t) == 2
                U2 = [U2, e];
            else
                U3 = [U3, e];
            end
        end

        % sizes must remain 5
        if numel(U1)~=5 || numel(U2)~=5 || numel(U3)~=5
            continue;
        end

        % sort within blocks
        U1 = sort(U1); U2 = sort(U2); U3 = sort(U3);

        % check each block is still a spanning tree
        if is_spanning_tree_edges(U1,n,edge_verts) && ...
           is_spanning_tree_edges(U2,n,edge_verts) && ...
           is_spanning_tree_edges(U3,n,edge_verts)

            % ---- KEY DIFFERENCE: canonicalize to ORDERED block order ----
            w = canonicalize_blocks_3_K6(U1, U2, U3);
            found = true;
            break;
        end
    end

    if ~found
        error('No valid ordered triangle-move w found for this (v,x,y,z).');
    end
end

% ===== canonicalize blocks so i1 < j1 < k1 =====
function row = canonicalize_blocks_3_K6(A, B, C)
    % A,B,C are already sorted internally
    firsts = [A(1), B(1), C(1)];
    [~, idx] = sort(firsts);           % reorder blocks by first element
    blocks = {A, B, C};
    row = [blocks{idx(1)}, blocks{idx(2)}, blocks{idx(3)}];
end

% ===== helper: lexicographic edge label for K_n =====
function lab = edge_label(u, v, n)
    if u > v
        tmp=u; u=v; v=tmp;
    end
    lab = sum(n - (1:(u-1))) + (v - u);
end

% ===== edge table for K_n in lex order =====
function edge_verts = edge_vertex_table(n)
    m = n*(n-1)/2;
    edge_verts = zeros(m,2);
    idx = 0;
    for a = 1:n
        for b = a+1:n
            idx = idx + 1;
            edge_verts(idx,:) = [a b];
        end
    end
end

% ===== check if a list of n-1 edges forms a spanning tree =====
function ok = is_spanning_tree_edges(edgeIdx, n, edge_verts)
    if numel(edgeIdx) ~= (n-1)
        ok = false; return;
    end
    s = zeros(numel(edgeIdx),1);
    t = zeros(numel(edgeIdx),1);
    for k = 1:numel(edgeIdx)
        uv = edge_verts(edgeIdx(k),:);
        s(k) = uv(1);
        t(k) = uv(2);
    end
    G = graph(s,t,[],n);
    ok = (numedges(G) == n-1) && (numel(unique(conncomp(G))) == 1);
end

% ===== check whole partition =====
function tf = is_tree_partition_K6(v)
    n = 6;
    if ~isvector(v) || numel(v) ~= 15
        tf = false; return;
    end
    v = v(:).';
    if ~isequal(sort(v), 1:15)
        tf = false; return;
    end
    edge_verts = edge_vertex_table(n);
    T1 = v(1:5); T2 = v(6:10); T3 = v(11:15);
    tf = is_spanning_tree_edges(T1,n,edge_verts) && ...
         is_spanning_tree_edges(T2,n,edge_verts) && ...
         is_spanning_tree_edges(T3,n,edge_verts);
end


Next, we have the code that gives us the permutations, in essence the elements of $H_3$.

function pxyz = triangle_move_permutation_K6_ordered(tree_partitions_file, x, y, z)
%TRIANGLE_MOVE_PERMUTATION_K6_ORDERED
% Permutation induced by a triangle move on ORDERED tree partitions of K6.

    % ---- validate x,y,z ----
    if ~(isscalar(x) && isscalar(y) && isscalar(z) && x < y && y < z)
        error('Need scalars x<y<z.');
    end
    if any([x y z] < 1) || any([x y z] > 6)
        error('Need 1 <= x < y < z <= 6.');
    end

    % ---- read partitions ----
    if ~ischar(tree_partitions_file) && ~isstring(tree_partitions_file)
        error('tree_partitions_file must be a filename (char or string).');
    end
    tree_partitions = readmatrix(tree_partitions_file);

    if size(tree_partitions,2) ~= 15
        error('File must have 15 columns per row (K6 ORDERED tree partitions).');
    end

    N = size(tree_partitions,1);

    % ---- build fast lookup: row -> index ----
    keys = strings(N,1);
    for i = 1:N
        keys(i) = row_key(tree_partitions(i,:));
    end
    dict = containers.Map(keys, 1:N);

    % ---- compute permutation ----
    pxyz = zeros(1, N);
    for k = 1:N
        v = tree_partitions(k,:);
        w = triangle_move_tree_partition_K6_ordered(v, x, y, z);

        kw = row_key(w);
        if ~isKey(dict, kw)
            error('Image partition not found in ORDERED file for row %d.', k);
        end
        pxyz(k) = dict(kw);
    end
end

function w = triangle_move_tree_partition_K6_ordered(v, x, y, z)
%TRIANGLE_MOVE_TREE_PARTITION_K6_ORDERED
% Triangle move, then canonicalize output to ORDERED convention:
% blocks sorted internally AND ordered by first entries i1<j1<k1.

    % ---- validate v ----
    if ~isvector(v) || numel(v) ~= 15
        error('v must be a 1x15 (or 15x1) vector.');
    end
    v = v(:).';

    if ~isequal(sort(v), 1:15)
        error('v must contain each label 1..15 exactly once.');
    end

    % ---- validate x,y,z ----
    if ~(isscalar(x) && isscalar(y) && isscalar(z) && x<y && y<z)
        error('Need scalars x<y<z.');
    end
    if any([x y z] < 1) || any([x y z] > 6)
        error('For K6, vertices must satisfy 1 <= x<y<z <= 6.');
    end

    if ~is_tree_partition_K6(v)
        error('Input v is not a tree partition of K6.');
    end

    n = 6;
    triEdges = [edge_label(x,y,n), edge_label(x,z,n), edge_label(y,z,n)];

    % split v into 3 blocks
    T1 = v(1:5);
    T2 = v(6:10);
    T3 = v(11:15);

    % mem(t) = which block contains triEdges(t)
    mem = zeros(1,3);
    for t = 1:3
        e = triEdges(t);
        if ismember(e,T1)
            mem(t)=1;
        elseif ismember(e,T2)
            mem(t)=2;
        else
            mem(t)=3;
        end
    end

    permsIdx = perms(1:3);          % 6 permutations
    edge_verts = edge_vertex_table(n);

    found = false;
    w = [];

    for r = 1:size(permsIdx,1)
        newMem = mem(permsIdx(r,:));

        % Must differ on at least TWO of the three triangle edges
        if sum(newMem ~= mem) < 2
            continue;
        end

        % remove tri edges from all blocks
        U1 = T1; U2 = T2; U3 = T3;
        for t = 1:3
            e = triEdges(t);
            U1(U1==e) = [];
            U2(U2==e) = [];
            U3(U3==e) = [];
        end

        % add them back according to newMem
        for t = 1:3
            e = triEdges(t);
            if newMem(t) == 1
                U1 = [U1, e];
            elseif newMem(t) == 2
                U2 = [U2, e];
            else
                U3 = [U3, e];
            end
        end

        if numel(U1)~=5 || numel(U2)~=5 || numel(U3)~=5
            continue;
        end

        % sort within blocks
        U1 = sort(U1); U2 = sort(U2); U3 = sort(U3);

        % check spanning trees
        if is_spanning_tree_edges(U1,n,edge_verts) && ...
           is_spanning_tree_edges(U2,n,edge_verts) && ...
           is_spanning_tree_edges(U3,n,edge_verts)

            % canonicalize block order so first entries are increasing
            w = canonicalize_blocks_3_K6(U1, U2, U3);
            found = true;
            break;
        end
    end

    if ~found
        error('No valid ordered triangle-move w found for this (v,x,y,z).');
    end
end

function row = canonicalize_blocks_3_K6(A, B, C)
% blocks are already internally sorted
    firsts = [A(1), B(1), C(1)];
    [~, idx] = sort(firsts);
    blocks = {A, B, C};
    row = [blocks{idx(1)}, blocks{idx(2)}, blocks{idx(3)}];
end

function lab = edge_label(u, v, n)
% Lexicographic edge label for K_n (u<v)
    if u > v
        tmp=u; u=v; v=tmp;
    end
    lab = sum(n - (1:(u-1))) + (v - u);
end

function edge_verts = edge_vertex_table(n)
% Edge list in lex order, rows are [u v]
    m = n*(n-1)/2;
    edge_verts = zeros(m,2);
    idx = 0;
    for a = 1:n
        for b = a+1:n
            idx = idx + 1;
            edge_verts(idx,:) = [a b];
        end
    end
end

function ok = is_spanning_tree_edges(edgeIdx, n, edge_verts)
% Check if edgeIdx (length n-1) forms a connected tree on n vertices.
    if numel(edgeIdx) ~= (n-1)
        ok = false; return;
    end
    s = zeros(numel(edgeIdx),1);
    t = zeros(numel(edgeIdx),1);
    for k = 1:numel(edgeIdx)
        uv = edge_verts(edgeIdx(k),:);
        s(k) = uv(1);
        t(k) = uv(2);
    end
    G = graph(s,t,[],n);
    ok = (numedges(G) == n-1) && (numel(unique(conncomp(G))) == 1);
end

function tf = is_tree_partition_K6(v)
% Check v is a partition of {1,...,15} and each block is a spanning tree.
    n = 6;
    if ~isvector(v) || numel(v) ~= 15
        tf = false; return;
    end
    v = v(:).';
    if ~isequal(sort(v), 1:15)
        tf = false; return;
    end
    edge_verts = edge_vertex_table(n);
    T1 = v(1:5); T2 = v(6:10); T3 = v(11:15);
    tf = is_spanning_tree_edges(T1,n,edge_verts) && ...
         is_spanning_tree_edges(T2,n,edge_verts) && ...
         is_spanning_tree_edges(T3,n,edge_verts);
end

function k = row_key(row)
% Unique key for matching rows in a hash map
    row = row(:).';
    k = strjoin(string(row), ',');
end

\noindent In this section, we give the MATLAB code for $G_3$. 

\noindent The following code gives you the ordered cycle-free partitions of $K_6$

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

Next, we have the code that describes the group action as seen in Lemma \ref{lemma 2}

function w = triangle_move_tree_partition_K6(v, x, y, z)
%TRIANGLE_MOVE_TREE_PARTITION_K6
% Given a tree partition v of K6 and vertices x<y<z in {1,...,6},
% produce a new tree partition w that differs from v only on the edges
% xy, xz, yz, and differs on at least TWO of those three edges.
%
% Representation:
%   v = [T1(1..5)  T2(1..5)  T3(1..5)]  (ordered blocks)
% where each Ti is a spanning tree on vertices 1..6.
%
% Edge labels in K6 (lex order):
% 1=(1,2), 2=(1,3), 3=(1,4), 4=(1,5), 5=(1,6),
% 6=(2,3), 7=(2,4), 8=(2,5), 9=(2,6),
% 10=(3,4),11=(3,5),12=(3,6),
% 13=(4,5),14=(4,6),
% 15=(5,6).

    % ---- validate inputs ----
    if ~isvector(v) || numel(v) ~= 15
        error('v must be a 1x15 (or 15x1) vector.');
    end
    v = v(:).'; % row

    if ~(isscalar(x) && isscalar(y) && isscalar(z))
        error('x, y, z must be scalars.');
    end
    if ~(x < y && y < z)
        error('Need x<y<z.');
    end
    if any([x y z] < 1) || any([x y z] > 6)
        error('For K6, vertices must satisfy 1 <= x<y<z <= 6.');
    end

    if ~isequal(sort(v), 1:15)
        error('v must contain each label 1..15 exactly once.');
    end

    if ~is_tree_partition_K6(v)
        error('Input v is not a tree partition of K6.');
    end

    % ---- triangle edges ----
    n = 6;
    triEdges = [edge_label(x,y,n), edge_label(x,z,n), edge_label(y,z,n)];

    % ---- split v into 3 trees (ordered blocks) ----
    T1 = v(1:5);
    T2 = v(6:10);
    T3 = v(11:15);

    % mem(t) in {1,2,3} tells which tree currently contains triEdges(t)
    mem = zeros(1,3);
    for t = 1:3
        e = triEdges(t);
        if ismember(e,T1), mem(t)=1;
        elseif ismember(e,T2), mem(t)=2;
        else, mem(t)=3;
        end
    end

    % We only change the membership of these 3 edges, keeping all others fixed.
    % To keep block sizes 5,5,5, each tree must receive the SAME COUNT of triangle edges as before.
    % So newMem must be a permutation of mem (same multiset of labels).
    permsIdx = perms(1:3);  % 6 permutations of positions (small, fine)
    found = false;
    w = [];

    edge_verts = edge_vertex_table(n); % for checking spanning trees quickly

    for r = 1:size(permsIdx,1)
        newMem = mem(permsIdx(r,:));  % permute assignments among the 3 edges

        % Must differ on at least TWO of the three triangle edges
        if sum(newMem ~= mem) < 2
            continue;
        end

        % Build new trees by removing all three triangle edges, then adding back according to newMem
        U1 = T1; U2 = T2; U3 = T3;

        % remove tri edges from wherever they are
        for t = 1:3
            e = triEdges(t);
            U1(U1==e) = [];
            U2(U2==e) = [];
            U3(U3==e) = [];
        end

        % add them back according to newMem
        for t = 1:3
            e = triEdges(t);
            if newMem(t) == 1
                U1 = [U1, e];
            elseif newMem(t) == 2
                U2 = [U2, e];
            else
                U3 = [U3, e];
            end
        end

        % sizes must remain 5
        if numel(U1)~=5 || numel(U2)~=5 || numel(U3)~=5
            continue;
        end

        % optional: sort within blocks (keeps your format stable)
        U1 = sort(U1); U2 = sort(U2); U3 = sort(U3);
        candidate = [U1 U2 U3];

        % check each block is still a spanning tree
        if is_spanning_tree_edges(U1,n,edge_verts) && ...
           is_spanning_tree_edges(U2,n,edge_verts) && ...
           is_spanning_tree_edges(U3,n,edge_verts)
            w = candidate;
            found = true;
            break;
        end
    end

    if ~found
        error('No valid triangle-move w found for this (v,x,y,z) under the "only triangle edges change" rule.');
    end
end

% ===== helper: lexicographic edge label for K_n =====
function lab = edge_label(u, v, n)
% u < v in {1,...,n}
    if u > v
        tmp=u; u=v; v=tmp;
    end
    lab = sum(n - (1:(u-1))) + (v - u);
end

% ===== edge table for K_n in lex order =====
function edge_verts = edge_vertex_table(n)
    m = n*(n-1)/2;
    edge_verts = zeros(m,2);
    idx = 0;
    for a = 1:n
        for b = a+1:n
            idx = idx + 1;
            edge_verts(idx,:) = [a b];
        end
    end
end

% ===== check if a list of n-1 edges forms a spanning tree =====
function ok = is_spanning_tree_edges(edgeIdx, n, edge_verts)
    if numel(edgeIdx) ~= (n-1)
        ok = false; return;
    end
    s = zeros(numel(edgeIdx),1);
    t = zeros(numel(edgeIdx),1);
    for k = 1:numel(edgeIdx)
        uv = edge_verts(edgeIdx(k),:);
        s(k) = uv(1);
        t(k) = uv(2);
    end
    G = graph(s,t,[],n);
    ok = (numedges(G) == n-1) && (numel(unique(conncomp(G))) == 1);
end

% ===== check whole partition =====
function tf = is_tree_partition_K6(v)
    n = 6;
    if ~isvector(v) || numel(v) ~= 15
        tf = false; return;
    end
    v = v(:).';
    if ~isequal(sort(v), 1:15)
        tf = false; return;
    end
    edge_verts = edge_vertex_table(n);
    T1 = v(1:5); T2 = v(6:10); T3 = v(11:15);
    tf = is_spanning_tree_edges(T1,n,edge_verts) && ...
         is_spanning_tree_edges(T2,n,edge_verts) && ...
         is_spanning_tree_edges(T3,n,edge_verts);
end

Next, we have the code that gives us the permutations, in essence the elements of $G_3$.
function pxyz = triangle_move_permutation_K6(tree_partitions_file, x, y, z)
%TRIANGLE_MOVE_PERMUTATION_K6  Permutation induced by a triangle move on K6.
%
% Inputs:
%   tree_partitions_file : txt file with K6 tree partitions (N x 15).
%                          Each row: [T1(1..5) T2(1..5) T3(1..5)]
%   x,y,z                : integers with 1 <= x < y < z <= 6
%
% Output:
%   pxyz : 1xN permutation vector. pxyz(k) is the row index of the image
%          of row k under triangle_move_tree_partition_K6(v,x,y,z).

    if nargin < 4
        error('Usage: pxyz = triangle_move_permutation_K6(tree_partitions_file, x, y, z)');
    end

    % ---- validate x,y,z ----
    if ~(isscalar(x) && isscalar(y) && isscalar(z) && x < y && y < z)
        error('Need scalars x<y<z.');
    end
    if any([x y z] < 1) || any([x y z] > 6)
        error('Need 1 <= x < y < z <= 6.');
    end

    % ---- read partitions ----
    if ~ischar(tree_partitions_file) && ~isstring(tree_partitions_file)
        error('tree_partitions_file must be a filename (char or string).');
    end

    tree_partitions = readmatrix(tree_partitions_file);

    if size(tree_partitions,2) ~= 15
        error('File must have 15 columns per row (K6 ordered tree partitions).');
    end

    N = size(tree_partitions,1);

    % ---- build fast lookup: row -> index ----
    keys = strings(N,1);
    for i = 1:N
        keys(i) = row_key(tree_partitions(i,:));
    end
    dict = containers.Map(keys, 1:N);

    % ---- compute permutation ----
    pxyz = zeros(1, N);
    for k = 1:N
        v = tree_partitions(k,:);
        w = triangle_move_tree_partition_K6(v, x, y, z);

        kw = row_key(w);
        if ~isKey(dict, kw)
            error('Image partition not found in file for row %d.', k);
        end
        pxyz(k) = dict(kw);
    end
end

function w = triangle_move_tree_partition_K6(v, x, y, z)
%TRIANGLE_MOVE_TREE_PARTITION_K6
% Changes ONLY the triangle edges xy,xz,yz, by reassigning them among the 3 blocks,
% keeping each block size = 5, and requiring the move changes >= 2 of the 3 triangle edges.

    % ---- validate v ----
    if ~isvector(v) || numel(v) ~= 15
        error('v must be a 1x15 (or 15x1) vector.');
    end
    v = v(:).';

    if ~isequal(sort(v), 1:15)
        error('v must contain each label 1..15 exactly once.');
    end

    % ---- validate x,y,z ----
    if ~(isscalar(x) && isscalar(y) && isscalar(z) && x<y && y<z)
        error('Need scalars x<y<z.');
    end
    if any([x y z] < 1) || any([x y z] > 6)
        error('For K6, vertices must satisfy 1 <= x<y<z <= 6.');
    end

    if ~is_tree_partition_K6(v)
        error('Input v is not a tree partition of K6.');
    end

    n = 6;
    triEdges = [edge_label(x,y,n), edge_label(x,z,n), edge_label(y,z,n)];

    % Split into 3 ordered blocks
    T1 = v(1:5);
    T2 = v(6:10);
    T3 = v(11:15);

    % mem(t) = which block contains triEdges(t)
    mem = zeros(1,3);
    for t = 1:3
        e = triEdges(t);
        if ismember(e,T1)
            mem(t)=1;
        elseif ismember(e,T2)
            mem(t)=2;
        else
            mem(t)=3;
        end
    end

    % Permute the assignments among the 3 triangle edges
    permsIdx = perms(1:3); % 6 permutations of positions
    edge_verts = edge_vertex_table(n);

    found = false;
    w = [];

    for r = 1:size(permsIdx,1)
        newMem = mem(permsIdx(r,:));

        % Must differ on at least TWO of the three triangle edges
        if sum(newMem ~= mem) < 2
            continue;
        end

        % Remove all triangle edges first
        U1 = T1; U2 = T2; U3 = T3;
        for t = 1:3
            e = triEdges(t);
            U1(U1==e) = [];
            U2(U2==e) = [];
            U3(U3==e) = [];
        end

        % Add them back according to newMem
        for t = 1:3
            e = triEdges(t);
            if newMem(t) == 1
                U1 = [U1, e];
            elseif newMem(t) == 2
                U2 = [U2, e];
            else
                U3 = [U3, e];
            end
        end

        if numel(U1)~=5 || numel(U2)~=5 || numel(U3)~=5
            continue;
        end

        % sort within blocks for stable matching to file rows
        U1 = sort(U1); U2 = sort(U2); U3 = sort(U3);
        candidate = [U1 U2 U3];

        if is_spanning_tree_edges(U1,n,edge_verts) && ...
           is_spanning_tree_edges(U2,n,edge_verts) && ...
           is_spanning_tree_edges(U3,n,edge_verts)
            w = candidate;
            found = true;
            break;
        end
    end

    if ~found
        error('No valid triangle-move w found for this (v,x,y,z).');
    end
end

function lab = edge_label(u, v, n)
% Lexicographic edge label for K_n (u<v).
    if u > v
        tmp=u; u=v; v=tmp;
    end
    lab = sum(n - (1:(u-1))) + (v - u);
end

function edge_verts = edge_vertex_table(n)
% Edge list in lex order, rows are [u v]
    m = n*(n-1)/2;
    edge_verts = zeros(m,2);
    idx = 0;
    for a = 1:n
        for b = a+1:n
            idx = idx + 1;
            edge_verts(idx,:) = [a b];
        end
    end
end

function ok = is_spanning_tree_edges(edgeIdx, n, edge_verts)
% Check if edgeIdx (length n-1) forms a connected acyclic graph on n vertices.
    if numel(edgeIdx) ~= (n-1)
        ok = false; return;
    end
    s = zeros(numel(edgeIdx),1);
    t = zeros(numel(edgeIdx),1);
    for k = 1:numel(edgeIdx)
        uv = edge_verts(edgeIdx(k),:);
        s(k) = uv(1);
        t(k) = uv(2);
    end
    G = graph(s,t,[],n);
    ok = (numedges(G) == n-1) && (numel(unique(conncomp(G))) == 1);
end

function tf = is_tree_partition_K6(v)
% Check each of the 3 blocks is a spanning tree and v is a partition of 1..15.
    n = 6;
    if ~isvector(v) || numel(v) ~= 15
        tf = false; return;
    end
    v = v(:).';
    if ~isequal(sort(v), 1:15)
        tf = false; return;
    end
    edge_verts = edge_vertex_table(n);
    T1 = v(1:5); T2 = v(6:10); T3 = v(11:15);
    tf = is_spanning_tree_edges(T1,n,edge_verts) && ...
         is_spanning_tree_edges(T2,n,edge_verts) && ...
         is_spanning_tree_edges(T3,n,edge_verts);
end

function k = row_key(row)
% Convert a numeric row vector into a unique string key for hashing.
    row = row(:).';
    k = strjoin(string(row), ',');
end

The following is the code that compares the orders of elements in $H_3$ and $G_3$.
function results = compare_tau_theta_generated_orders_verbose()
%COMPARE_TAU_THETA_GENERATED_ORDERS_VERBOSE
% Like your code, but when a mismatch occurs it prints:
%   - ALL order-keys where counts differ
%   - (optionally) a witness subset of indices for each differing key

    % -----------------------------
    % USER SETTINGS
    % -----------------------------
    KMAX = 20;                 % test k=2..KMAX
    PRINT_ALL_BAD_KEYS = true; % print every differing order-key at mismatch k
    MAX_WITNESSES = 25;        % how many witness subsets to print (set Inf for all)
                               % (printing all may be large)
    FIND_WITNESS_FOR_EACH_KEY = true;

    % >>>>> EDIT THESE LISTS <<<<<
    tau_files = {
        'k6permutations123.txt'
        'k6permutations124.txt'
        'k6permutations125.txt'
        'k6permutations126.txt'
        'k6permutations134.txt'
        'k6permutations135.txt'
        'k6permutations136.txt'
        'k6permutations145.txt'
        'k6permutations146.txt'
        'k6permutations156.txt'
        'k6permutations234.txt'
        'k6permutations235.txt'
        'k6permutations236.txt'
        'k6permutations245.txt'
        'k6permutations246.txt'
        'k6permutations256.txt'
        'k6permutations345.txt'
        'k6permutations346.txt'
        'k6permutations356.txt'
        'k6permutations456.txt'
    };

    theta_files = {
        'k6permutationsordered123.txt'
        'k6permutationsordered124.txt'
        'k6permutationsordered125.txt'
        'k6permutationsordered126.txt'
        'k6permutationsordered134.txt'
        'k6permutationsordered135.txt'
        'k6permutationsordered136.txt'
        'k6permutationsordered145.txt'
        'k6permutationsordered146.txt'
        'k6permutationsordered156.txt'
        'k6permutationsordered234.txt'
        'k6permutationsordered235.txt'
        'k6permutationsordered236.txt'
        'k6permutationsordered245.txt'
        'k6permutationsordered246.txt'
        'k6permutationsordered256.txt'
        'k6permutationsordered345.txt'
        'k6permutationsordered346.txt'
        'k6permutationsordered356.txt'
        'k6permutationsordered456.txt'
    };

    if numel(tau_files) ~= 20 || numel(theta_files) ~= 20
        error('Need exactly 20 tau files and 20 theta files.');
    end

    % -----------------------------
    % LOAD GENERATORS
    % -----------------------------
    fprintf('Loading tau generators...\n');
    tau = load_generators(tau_files);

    fprintf('Loading theta generators...\n');
    theta = load_generators(theta_files);

    nTau = numel(tau{1});
    nTheta = numel(theta{1});
    fprintf('tau in S_%d, theta in S_%d\n', nTau, nTheta);

    for i = 1:20
        check_is_permutation(tau{i});
        check_is_permutation(theta{i});
    end

    % -----------------------------
    % MAIN EXPERIMENT
    % -----------------------------
    results = struct();
    results.KMAX = KMAX;
    results.nTau = nTau;
    results.nTheta = nTheta;
    results.perK = cell(1, KMAX);
    results.firstMismatchK = NaN;
    results.badKeys = {};
    results.badKeyCounts = [];
    results.witnesses = {};

    gens = 1:20;

    for k = 2:KMAX
        fprintf('\n==== k = %d ====\n', k);

        combos = nchoosek(gens, k);
        numC = size(combos, 1);
        fprintf('Number of combinations: %d\n', numC);

        tauHist = containers.Map('KeyType','char','ValueType','double');
        thetaHist = containers.Map('KeyType','char','ValueType','double');

        for r = 1:numC
            idx = combos(r,:);

            Pt = product_of_generators(tau, idx);
            Pq = product_of_generators(theta, idx);

            keyTau = order_prime_factor_key(Pt);
            keyTheta = order_prime_factor_key(Pq);

            tauHist = map_increment(tauHist, keyTau);
            thetaHist = map_increment(thetaHist, keyTheta);

            if mod(r, 2000) == 0
                fprintf('  processed %d / %d\n', r, numC);
            end
        end

        same = maps_equal(tauHist, thetaHist);

        results.perK{k} = struct( ...
            'k', k, ...
            'numCombos', numC, ...
            'sameMultiset', same, ...
            'tauDistinctOrders', tauHist.Count, ...
            'thetaDistinctOrders', thetaHist.Count ...
        );

        fprintf('Distinct order-types: tau=%d, theta=%d\n', tauHist.Count, thetaHist.Count);
        fprintf('Multiset match at k=%d?  %s\n', k, string(same));

        if ~same
            fprintf('\n*** MISMATCH FOUND at k=%d ***\n', k);

            % ---- compute ALL bad keys ----
            [badKeys, badCounts] = all_bad_keys(tauHist, thetaHist);

            results.firstMismatchK = k;
            results.badKeys = badKeys;
            results.badKeyCounts = badCounts;

            fprintf('Number of differing order-keys: %d\n', numel(badKeys));

            if PRINT_ALL_BAD_KEYS
                fprintf('\nList of differing order-keys (key : tauCount vs thetaCount)\n');
                for i = 1:numel(badKeys)
                    fprintf('%s : %d vs %d\n', badKeys{i}, badCounts(i,1), badCounts(i,2));
                end
            end

            % ---- find witnesses (optional) ----
            if FIND_WITNESS_FOR_EACH_KEY
                fprintf('\nFinding witness subsets (up to %d shown)...\n', MAX_WITNESSES);

                witnesses = cell(numel(badKeys),1);
                shown = 0;

                for i = 1:numel(badKeys)
                    key = badKeys{i};

                    w = find_witness_for_key(tau, theta, combos, key);
                    witnesses{i} = w;

                    if ~isempty(w)
                        shown = shown + 1;
                        fprintf('\nWitness for key %s\n', key);
                        fprintf('subset indices: ');
                        disp(w.idx);
                        fprintf('tau key   = %s\n', w.tauKey);
                        fprintf('theta key = %s\n', w.thetaKey);

                        if shown >= MAX_WITNESSES
                            fprintf('\nStopped printing witnesses after %d.\n', MAX_WITNESSES);
                            break;
                        end
                    end
                end

                results.witnesses = witnesses;
            end

            return;
        end
    end

    fprintf('\nAll k=2..%d matched under this test.\n', KMAX);
end

% ============================================================
% BAD KEY COMPUTATION
% ============================================================
function [badKeys, badCounts] = all_bad_keys(A, B)
% Return all keys where counts differ. badCounts(i,:) = [countA countB]

    keysA = keys(A);
    keysB = keys(B);

    % union of keys
    allKeys = unique([string(keysA(:)); string(keysB(:))]);

    badKeys = {};
    badCounts = zeros(0,2);

    for i = 1:numel(allKeys)
        k = char(allKeys(i));

        ca = 0; cb = 0;
        if isKey(A,k), ca = A(k); end
        if isKey(B,k), cb = B(k); end

        if ca ~= cb
            badKeys{end+1,1} = k; %#ok<AGROW>
            badCounts(end+1,:) = [ca cb]; %#ok<AGROW>
        end
    end
end

function w = find_witness_for_key(tau, theta, combos, key)
% Find one subset idx where tau gives 'key' but theta doesn't,
% OR theta gives 'key' but tau doesn't.

    w = struct('idx',[],'tauKey','','thetaKey','');

    for r = 1:size(combos,1)
        idx = combos(r,:);

        Pt = product_of_generators(tau, idx);
        kt = order_prime_factor_key(Pt);

        Pq = product_of_generators(theta, idx);
        kq = order_prime_factor_key(Pq);

        isTau = strcmp(kt, key);
        isTheta = strcmp(kq, key);

        if isTau ~= isTheta
            w.idx = idx;
            w.tauKey = kt;
            w.thetaKey = kq;
            return;
        end
    end
end

% ============================================================
% Helper functions
% ============================================================

function gens = load_generators(files)
    gens = cell(1, numel(files));
    for i = 1:numel(files)
        v = readmatrix(files{i});
        gens{i} = v(:).';
    end
end

function check_is_permutation(p)
    n = numel(p);
    if ~isequal(sort(p), 1:n)
        error('A generator is not a permutation of 1..n.');
    end
end

function w = compose(u, v)
    w = u(v);
end

function P = product_of_generators(gensCell, idx)
    P = gensCell{idx(1)};
    for t = 2:numel(idx)
        P = compose(P, gensCell{idx(t)});
    end
end

function key = order_prime_factor_key(p)
    L = permutation_cycle_lengths(p);
    if isempty(L)
        key = '1';
        return;
    end
    key = lcm_prime_factor_key_from_lengths(L);
end

function cycle_lengths = permutation_cycle_lengths(p)
    n = numel(p);
    visited = false(1,n);
    cycle_lengths = zeros(1,n);
    ccount = 0;

    for i = 1:n
        if ~visited(i)
            len = 0;
            j = i;
            while ~visited(j)
                visited(j) = true;
                j = p(j);
                len = len + 1;
            end
            if len > 1
                ccount = ccount + 1;
                cycle_lengths(ccount) = len;
            end
        end
    end
    cycle_lengths = cycle_lengths(1:ccount);
end

function key = lcm_prime_factor_key_from_lengths(L)
    M = max(L);
    plist = primes(M);

    parts = strings(0,1);
    for p = plist
        e = 0;
        pe = p;
        while pe <= M
            if any(mod(L, pe) == 0)
                e = e + 1;
                if pe > floor(M / p)
                    break;
                end
                pe = pe * p;
            else
                break;
            end
        end

        if e > 0
            if e == 1
                parts(end+1) = string(p); %#ok<AGROW>
            else
                parts(end+1) = string(p) + "^" + string(e); %#ok<AGROW>
            end
        end
    end

    if isempty(parts)
        key = '1';
    else
        key = char(strjoin(parts, '*'));
    end
end

function m = map_increment(m, key)
    if isKey(m, key)
        m(key) = m(key) + 1;
    else
        m(key) = 1;
    end
end

function tf = maps_equal(A, B)
    if A.Count ~= B.Count
        tf = false; return;
    end
    kA = keys(A);
    for i = 1:numel(kA)
        k = kA{i};
        if ~isKey(B, k)
            tf = false; return;
        end
        if A(k) ~= B(k)
            tf = false; return;
        end
    end
    tf = true;
end



