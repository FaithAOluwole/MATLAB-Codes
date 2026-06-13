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