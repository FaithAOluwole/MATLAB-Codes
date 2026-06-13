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