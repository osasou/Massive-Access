classdef Decoder
%Class responsible for decoding chirps - contains all functons pertaining to decoding
    properties
        Y       % Y{p} is a 2^m x 2^p matrix of measurements for patch p
        r       % 0, 1 or 2; 2^r patches
        re      % logical: false=complex chirps, true=real chirps
        m       % size of chirp code (2^m is length of codewords in each slot) recommend m = 6, 7 or 8
        p       % 2^p slots, require p<=(m-r)(m-r+3)/2-1 for complex
                %                    p<=(m-r)(m-r+1)/2-1 for real
        K       % number of messages
        params  % various; see below
        patches % number of patches
    end

    methods
        function self = Decoder(Y,re,m,p,K)
            addpath('utils');
            self.Y=Y;
            self.r=0;
            self.re=re;
            self.m=m;
            self.p=p;
            self.K=K;
            self.patches=2^self.r;
            
            %default parameter values
            self.params.alpha = 0.1; %accept components if coefficient is within alpha of 1
            
            %expected sparsity for a given slot
            self.params.tree_order = 1; %3; %number of candidates per row in findPb tree search
            
        end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [output_bits, h_hat] = chirrup_decode(self)

        % chirrup_decode  Performs CHIRRUP decoding
        %
        % output_bits  B x K matrix of the K B-bit messages
        % timing       running time in seconds
        %
        % AJT (12/9/18)

            global outer_recov
            warning('off','MATLAB:rankDeficientMatrix');

            for patch = 1:self.patches %patchは一つでやる

                    outer_recov = [];
                    Yp = self.Y{patch}; % 1patch分のY、つまり2^p slotsはある

                    %cycle through the slots repeatedly
                    for slot = 1:2^self.p %slotごとに計算
                        %run chirp reconstruction on given slot
                        recov = self.chirp_rec(Yp(:,slot),self.K,slot);                            
                    end

                    outer_recov = recov;
                    
                    %convert components back into bits
                    bit_matrix = self.find_bits(outer_recov);
                    
                    h_hat=[];
                    for i=1:length(outer_recov)
                        h_hat = [h_hat, outer_recov(i).c];
                    end
            end
                    output_bits=[];
                    output_bits = bit_matrix';
                    kmin = min(self.K,size(output_bits,2));
                    h_hat = h_hat(:,1:kmin);
                    output_bits = output_bits(:,1:kmin);
        end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        function recov = chirp_rec(self,y,nitlim,slot)

        % chirp_rec     Runs the chirp reconstruction algorithm,
        %               incorporating(含んでいる) components found in other slots.
        %
        % y        measurement vector for given slot
        % nitlim   iteration limit
        % slot     slot number
        %
        % recov    multidimensional struct storing (P,b,c) for each component
        %          e.g. for component i: recov(1).P, recov(1).b, recov(1).c
        %
        % Sina Jafarpour (14/4/08) and AJT (12/9/18)

            global outer_recov

            foo.P = [];
            foo.b = [];
            foo.c = [];
            foo.comps = [];
            recov(1) = foo;
            allRM = [];
            ncomp = 0;
            user_num = self.K;

            while (ncomp < user_num && norm(y)>1e-3)

                [Phat, bhat, y_1, y_2, b_1, beta_1] = self.findPb(y);

                %determine component from P,b
                RM = Encoder.gen_chirp(Phat,bhat);
                allRM = [allRM RM];

                %use all components & refine previous ests.
                h = 0.5 * (y_1 + conj((-1)^(b_1 + 0.5 * beta_1)) * y_2);
                y = y - RM*h;

                %a component has been detected and handled.
                %record component parameters
                ncomp = ncomp+1;
                recov(ncomp).P = Phat;
                recov(ncomp).b = bhat;
                recov(ncomp).c = h;
                recov(ncomp).comp = RM;

            end
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function bits = find_bits(self,res)

            bits = [];
            K = length(res);

            for k = 1:K

                Pk = res(k).P;
                bk = res(k).b;
                [i, j] = upper_indices(self.m, self.re);
                for l = 1:length(i)
                    Pvec(l) = Pk(i(l),j(l));
                end
                bits2 = [Pvec bk];
                bits(k,:) = [bits2];
            end
        end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [Pout, bout, y_1, y_2, b_1, beta_1] = findPb(self,y)

            alpha_all = [];
            b = [];
            beta = [];
            
            
            twoM = length(y);
            M = log2(twoM);
            twom = length(y);
            m = log2(twom);
            y_loop = y;
            s = twom;
            
            while m > 1
%                 disp(s)
                s = s - 0.5*s;

                y_even = [];
                y_odd = [];
                
                for i = 2:2:size(y_loop)
                    y_even = [y_even y_loop(i)];
                end
                
                for i = 1:2:size(y_loop)
                    y_odd = [y_odd y_loop(i)];
                end

                y_conj = y_even.*conj(y_odd);

                y_fht = self.fhtnat(y_conj');
                
                if m > 2
                    cand = zeros(2^(m - 1), m - 1);
                    A = [0:2^(m - 1) - 1];
                    bit = de2bi(A);
                    bits = flip(bit);
                    bits = flip(bit, 2);
                    cand(:,1:m - 1) = bits;
                else
                    cand = [0, 1].';
                end


                [foo, ind] = sort(-abs(y_fht));

                alpha_all = [alpha_all, cand(ind(1),:)];
                alpha = cand(ind(1),:);

                max = y_fht(ind(1));
                max_re = real(max);
                max_imag = imag(max);
                if max_imag == 0 && max_re == 0
                    b = [b, 0];
                    beta = [beta, 0];
                    b_now = 0;
                    beta_now = 0;
                    disp("zero");
                end
                if abs(max_imag) < abs(max_re)
                    
                    if max_re < 0
                        b = [b, 1];
                        beta = [beta, 0];
                        b_now = 1;
                        beta_now = 0;
                    end
                    
                    if max_re > 0
                        b = [b, 0];
                        beta = [beta, 0];
                        b_now = 0;
                        beta_now = 0;
                    end
                
                else
                    
                    if max_imag < 0
                        b = [b, 1];
                        beta = [beta, 1];
                        b_now = 1;
                        beta_now = 1;
                    end

                    if max_imag > 0
                        b = [b, 0];
                        beta = [beta, 1];
                        b_now = 0;
                        beta_now = 1;
                    end
                end

                v = (-1).^(b_now+beta_now/2+(alpha*cand'));

                y_loop = 0.5 .* (y_odd + conj(v) .* y_even);
                y_loop = y_loop.';
                m = m - 1;
            end
            
            y_1 = y_loop(1);
            y_2 = y_loop(2);
            y_zero = conj(y_loop(1))*y_loop(2);
            
            y_zero_re = real(y_zero);
            y_zero_imag = imag(y_zero);
            
            if y_zero_imag == 0 && y_zero_re == 0
                b = [b, 0];
                beta = [beta, 0];
                disp("zero");
            end
            
            if abs(y_zero_imag) < abs(y_zero_re)
                if y_zero_re < 0
                    b = [b, 1];
                    beta = [beta, 0];
                end
                if y_zero_re > 0
                    b = [b, 0];
                    beta = [beta, 0];
                end
            else
                if y_zero_imag < 0
                    b = [b, 1];
                    beta = [beta, 1];
                end
                if y_zero_imag > 0
                    b = [b, 0];
                    beta = [beta, 1];
                end
            end
            
            m = log2(twom);
            alpha_all = fliplr(alpha_all); %% ここ変える
            b = fliplr(b);
            beta = fliplr(beta);
            P = [beta(1)];
            b_1 = b(1);
            beta_1 = beta(1);

            cnt = 0;
            pos = 1;
            for s = 2:m
                
                P_foo = zeros(s,s);
                P_foo(1:s-1, 1:s-1) = P;
                P_foo(1:1+cnt, s) = alpha_all(pos:pos+cnt);
                P_foo(s, 1:1+cnt) = alpha_all(pos:pos+cnt);
                P_foo(s,s) = beta(s);
                pos = pos + cnt + 1;
                cnt = cnt + 1;
                P = P_foo;
%                 disp(P_foo)
            end
            
            Pout = P;
            bout = b;

        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [Pout,bout,crec] = fPN_recursive(self,m,y,P,rec)

            global perm prout DGbasis;
            %defined prior to call of this recursive function

            M = length(P);

            %which basis matrices go with current row
            DGblen = size(DGbasis,3);
            if (self.re==0)
                csM1 = cumsum(M:-1:1);
            else
                csM1 = cumsum(M-1:-1:1);
            end
            DGto = csM1(m);
            DGto = min(DGto,DGblen);

            pm = perm(m);
            row = P(pm,:);
            if (self.re==0)
                fixed = perm(1:m-1);
                unfixed = perm(m:end);
            else
                fixed = perm(1:m);
                unfixed = perm(m+1:end);
            end

            %generate all candidates for this row,
            %given constraints by already-fixed rows
            %(= fixed element values in this row)
            nfree = 2^length(unfixed);

            nums = 0:nfree-1;
            cands = row(ones(nfree,1),:);
            cands(:,unfixed) = (dec2bin(nums)=='1');
            candix = cands * (2.^(M-1:-1:0)')  +1;  % note +1 for Matlab index

            %probe values of these candidates
            %when probed with vector having 1 in pm_th pos
            s1 = abs(prout(pm,candix));

            if m>1

                % use vector with ones in current pos plus previous pos
                pmprev = perm(m-1);
                evpch = zeros(M,1);  evpch([pm pmprev]) = 1;
                prout_pch = self.probe_ye(y,evpch);
                ixpch = rec(m-1).candix-1;  % note -1 coz it was a matlab index.
                ixvec = bitxor(candix-1,ixpch);
                s2 = abs(prout_pch(ixvec+1))';

            else
                s2 = zeros(size(s1));
                %pmprev = [];
            end

            % order candidates based on scores
            [foo,rix] = sort(-(s1+s2));    % score + parity score

            rmax = s1(rix);
            rmax_pch = s2(rix);

            % this is used in b computation.
            % Set up here to avoid doing in every loop
            a = (0:2^M-1)';
            aval = double( dec2bin(a,M)=='1' );  % binary expansion

            % process each candidate in order
            Pout = [];
            crec = [];
            bout = [];

            for v = 1:min(nfree,self.params.tree_order)

                vsc = rmax(v);
                vix = rix(v);
                psc = rmax_pch(v);

                newrow = cands(vix,:);
                Phat = P;
                Phat(pm,:) = newrow;
                Phat(:,pm) = Phat(pm,:);

                newrec.candix = candix(vix);
                newrec.vsc = vsc;
                newrec.psc = psc;

                recx = [rec; newrec];

                %have we filled those rows of matrix which are constrained by
                %the D-G basis?
                if DGto==(DGblen)

                    Psynth = Phat;

                    % "detection"
                    % compute Hadamard probe on the data for each row of P matrix.
                    % accumulate Hadamard power and perform Sequential test until
                    % decision on presence/absence of component is made.
                    hadpow = zeros(2^M,1);
                    detected = 0;
                    Mvec = 0:2^M-1;
                    for r = 1:M
                        % probe data with this row of matrix.
                        ev = 2^(r-1);
                        evb = abs(dec2bin(ev,M)=='1');
                        [hx,yp] = self.probe_ye(y,ev);
                        % permute wrt P row, so peak should be at index zero
                        prows = mod( sum(Psynth(M-r+1,:),1), 2);
                        % note M-r+1 should really be find(ev)
                        % but I can cheat 'coz I know location of only 1-value
                        rval = (2.^(M-1:-1:0)) * prows';
                        ixvec = bitxor(Mvec,rval);
                        hxs = hx(ixvec+1);
                        % add on this power
                        hadpow = hadpow + real(hxs.*conj(hxs));
                        foo = self.findoutliers(hadpow, 2.5);
                        if foo(1)>0 % was ~=0 but we want big outliers only
                            detected = 1;
                        end

                        if detected~=0
                            break;  % decision made so stop loop
                        end
                    end
                    if detected==1
                        % confirm detection and find b
                        % via HadXform
                        Phat = Psynth;
                        euse = zeros(M,1);
                        Pe = rem( Phat*euse ,2);
                        phi = Encoder.gen_chirp(Phat,Pe');
                        foo = conj(phi).*y.*( (-1).^( aval*Pe) );
                        Hfoo = self.fhtnat(foo);  % * scale 2^M
                        Hfoo = Hfoo * 2^M;  % restore scale
                        [qq,ww] = max(abs(Hfoo));
                        bhat = aval(ww,:);

                        % now test this b.
                        % power in max of Had xform
                        pmax = qq^2;
                        % power in rest
                        ptot = real(Hfoo'*Hfoo);
                        prest = ptot-pmax;
                        % power per bin
                        prestbin = prest/(2^M-1);
                        % would expect this to be power of bin sans component.
                        % and thus standard deviation of amplitude should be
                        sigamp = sqrt(prestbin);
                        confirmed = qq>3*sigamp;

                        %accept candidate only if b test is good.
                        if confirmed
                            Pout = Phat;
                            crec = recx;
                            bout = bhat;
                            break;  % terminate loop;  accept this candidate.
                        end
                    end
                else

                    % recursion for next row
                    [Pout_, bout_, crec_] = self.fPN_recursive(m+1,y,Phat,recx);

                    if ~isempty(Pout_)
                        % if we get here, return cand deemed OK
                        Pout = Pout_;
                        bout = bout_;
                        crec = crec_;
                        break;  % quit loop & return
                    end

                end % if DGto==DGblen
            end  % for v

        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [prb,yprod] = probe_ye(self,y,ein)

        % probe_ye  probe data vector with error vector
        %
        % given data vector  "" y ""  and "error value"   "" e "",
        % compute Hadamard transform ""prb""  of yprod = conj(y(a)) * y(a+e)


        % Sina Jafarpour 28/3/08

            Mpow = length(y);
            M = log2(Mpow);  % this should be an integer
            %isscalar(A):Aがスカラーなら1、その他のときは0
            %einはevp(単位行列)の上から一行ずつくる
            if isscalar(ein)
                e = ein;
            else
                e = (2.^(M-1:-1:0))*ein(:);
                %ein(:)はein()を列ベクトルに形状変更
                %2.^(M-1:-1:0)は [2^(M-1) 2^(M-2) 2^(M-3) ... 2^(0)]となる
                %つまり、[2^(M-1) 2^(M-2) 2^(M-3) ...2^(0)]・ein^Tのベクトルの内積をとっている
            end
            %{
            disp(ein(:))
            disp("2^  ")
            disp(2.^(M-1:-1:0))
            %}

            % index vector (starting at 0)
            a = (0:Mpow-1)';
            % a : [0,1,2,3,...,Mpow-1]^T 列ベクトル

            % a plus e values
            ape = bitxor(a,e);

            % y(a+e)
            yape = y(ape+1); % ← どんな計算しとるん？？？？？ 1×nの列ベクトルと1×nの列ベクトルの掛け算？
            % シフトしている！！！
            
            % product
            %conj(x):xの複素共役を求める
            yprod = conj(y).*yape;

            % use fast Walsh-Hadamard transform routine
            prb = self.fhtnat(yprod);

        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function x=fhtnat(self, data)
        %------------------------------------------------------
        %1D Natural(Hadamard)ordered Fast Hadamard Transform
        %------------------------------------------------------
        % Author: Gylson Thomas
        % e-mail: gylson_thomas@yahoo.com
        % Asst. Professor, Electrical and Electronics Engineering Dept.
        % MES College of Engineering Kuttippuram,
        % Kerala, India, February 2005.
        % copyright 2007.
        % This function implements the 1D natural(Hadamard)ordered Fast Hadamard Transform,
            N = pow2(floor(log2(length(data))));
            x = data(1:N);
            k1=N; k2=1; k3=N/2;
            for i1=1:log2(N)
                L1=1;
                for i2=1:k2
                    for i3=1:k3
                        i=i3+L1-1; j=i+k3;
                        temp1= x(i); temp2 = x(j);
                        x(i) = temp1 + temp2;
                        x(j) = temp1 - temp2;
                    end
                        L1=L1+k1;
                end
                    k1 = k1/2;  k2 = k2*2;  k3 = k3/2;
            end
            %inv(A) : Aの逆行列を示す
            x=inv(N)*x; %Delete this line for inverse transform
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [P,b] = bstr2Pb(self,M,bstr)
        % convert a bitstring (as an integer or as a logical vector) to P and b
        % * updated for Kerdock codes:  no longer zero diagonal.
        %
        % [P,b] = bstr2Pb(M,[]);                randomly assigned P and b
        % [P,b] = bstr2Pb(M,[1 0 1 0 ...])      take bits from vector elements
        % [P,b] = bstr2Pb(M,int);               take bits from int

        % began life as a private function of genRM.m
        % SJS 6/11/07
        % SJS 28/3/08 rehashed for Kerdock codes:
        %             diagonal no longer constrained to be all-zeros
        %             total # free elements thus M(M+3)/2

            if self.re == 0
                MM = ((M^2+3*M)/2); %    MMold = ((M^2+M)/2);
                ixmid = (M+1)*M/2;
            else
                MM = 0.5*M*(M+1);
                ixmid = (M-1)*M/2;
            end

            if isempty(bstr)
                %bstr = floor(rand*2^MM);  % precision errors!
                bstr = randn(1,MM)>0;
            end
            if isscalar(bstr)
                % could use dec2bin...
                bvec = false(1,MM);
                for q = 1:MM
                    bit = rem(bstr,2); % bitand(bstr,1);
                    bvec(end-q+1) = bit;
                    bstr = floor(bstr/2);
                end
            else
                bvec = bstr;
                % should check that length is M, but what to do if it ain't?
            end


            % fprintf('MM %i (prev %i)   ix %i (prev %i)   len %i\n',...
            %     MM,MMold,ixmid,ixmidold,length(bvec));
            P = false(M);
            if self.re == 0
                ix = find(tril(ones(M),0));  % indices of lower tri, including diag
            else
                ix = find(tril(ones(M),-1)); % doesnt include the diagonal
            end
            P(ix) = bvec(1:ixmid);
            P = P | P';  % make symmetric
            b = bvec(ixmid+1:end);

            end
        


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function out = findoutliers(self,dat,mth)

        % findoutliers  Find outliers in a single data vector.
        %
        % foo = findoutliers(dat,mth)
        %       dat must be a vector.
        %
        % foo = nonzero where dat element was an outlier
        %       where possible, sign of foo indicates "direction" of outlier
        %
        % mth = method:
        %       1   = 3 sigmas from the mean
        %       2   = quartile test, Mendenhall & Sincich quartiles
        %       2.5 = quartile test, Tukey quartiles
        %
        % see  http://mathforum.org/library/drmath/view/52720.html

            if nargin<2 | isempty(mth)
                mth = 2;
            end
            ndat = length(dat);

            switch floor(mth)
                case 1
                    % 3 standard deviations from mean
                    if max(size(dat))<=1
                        out = zeros(size(dat));
                        return
                    else
                        sd = max(std(dat),eps);
                    end

                    mu = mean(dat);
                    demu = abs(dat-mu)/sd;
                    out = demu>3;

                    % introduce sign in order to tell if greater or lesser.
                    out = out.*sign(dat-mu);

                case 2
                    [sdat,ix] = sort(dat);
                    [foo,undoix] = sort(ix);  % for unsorting
                    %pctile = (1:ndat)/ndat;

                    % compute quartiles
                    switch mth-floor(mth)
                        case 0
                            % Mendenhall and Sincich method
                            L = ceil(0.25*(ndat+1));
                            U = floor(0.75*(ndat+1));
                            LQ = sdat(L);
                            UQ = sdat(U);
                        case .5
                            % Tukey method
                            LQ = median(sdat(1:ceil(ndat/2)));
                            UQ = median(sdat(ceil((ndat+1)/2):ndat));
                    end

                    IQR = UQ-LQ;
                    out = (sdat>UQ+1.5*IQR) - (sdat<LQ-1.5*IQR);
                    out = out(undoix);
                    %fprintf('findoutliers:  LQ %g  UQ %g  IQR %g  no %g\n',LQ,UQ,IQR,sum(out~=0))

                otherwise
                    error(['no such method:  ' int2str(mth)])
            end
        end



    end
end