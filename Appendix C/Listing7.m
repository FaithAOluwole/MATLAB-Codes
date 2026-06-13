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