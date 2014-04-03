function [imf,ort] = memd_emd(varargin)
% standard EMD from [2,3] with some tweaks.
%
% SYNTAX
%   imf = memd_emd(s)
%   imf = memd_emd(s,...,'Option_name',Option_value,...)
%   imf = memd_emd(s,opt)
%   [imf,ort,order_der] = memd_emd(...)
%
% INPUTS
%   s : one-dimensional real signal
%   opt : a struct containing optional parameters. They are listed below
%           together with the default value.  <opt> can be given as an
%           already-created struct, or it can be given implicitly as a
%           series of fieldname-value pairs that are passed to the Matlab
%           'struct' function.
%       alpha       <0.05>  : Stopping criterion parameter [3].
%       localEMDfunc <[]>   : Handle to a function that accepts as
%                       arguments: the residual signal <r>, the envelope of
%                       maxima <envmax>, the envelope of minima <envmin>,
%                       the scalar parameter <alpha>, and a variable to
%                       hold any needed parameter values <localEMDparam>;
%                       and returns a vector of the same size containing
%                       the weighting function as described in [3] section
%                       3.3.  Since this technique can make it impossible
%                       to meet an independent stopping criterion, the
%                       stopping criterion is replaced by the condition
%                       that all elements of the weighting function are
%                       zero.
%       localEMDparam <{}> : a value containing any additional arguments
%                       to parameterize the <localEMDfunc>.  Do not use
%                       a cell array here because the call to 'struct' in
%                       'init' will do undesirable things with it.
%                       However, a struct will work fine for passing in
%                       multiple values of different types.
%       maxmodes    <8>     : Maximum number of IMFs; <0> means no maximum.
%       postprocess <empty> : Handle to a function that accepts as
%                       arguments (1) the newly extracted IMF <r> and (2)
%                       the structure of auxiliary values returned by its
%                       matching preprocess function.  If there is no
%                       preprocessing function, the second argument will be
%                       empty. The function is called immediately following
%                       the termination of sifting for each IMF. <r> is a
%                       row vector when passed in to the postprocess
%                       function, which must return a modified version as a
%                       row vector. This is where you would e.g. remove the
%                       masking signal added by <preprocess>.
%       preprocess  <empty> : Handle to a function that accepts as
%                       arguments the residual signal <r> and a value
%                       (which could be a struct, scalar, cell array, or
%                       whatever) containing any necessary parameters. It
%                       must return a modified version of <r> that is used
%                       in place of the original <r>, and a value
%                       containing any auxiliary data that will be needed
%                       by the matching 'postprocess'.  The preprocess
%                       function is executed immediately before beginning
%                       the sifting for each IMF.  <r> is a row vector when
%                       passed to the preprocess function, and the modified
%                       version returned must also be a row vector. This is
%                       where you would e.g. add a masking signal.
%       pre_params <empty>  : a value containing the parameters needed by
%                       the matching <preprocess>.
%       stop        <'f'>   : Kind of sifting stopping criterion : 'h' for
%                       the Huang criterion [2] and 'f' for Flandrin and
%                       Rilling's [3].
%       t                   : Time vector t
%
% OUTPUTS
%   imf     : (number of IMFs +1) X length(s) array which contains each
%               IMFs and a residue
%   ort     : orthogonality index defined in [3,1]
%
% REFERENCES
%   [1] T. Oberlin, S. Meignen and V. Perrier, An Alternative Formulation
%       for the Empirical Mode Decomposition. IEEE TRANSACTIONS ON SIGNAL
%       PROCESSING, VOL. 60, NO. 5, MAY 2012 
%   [2] N. Huang, Z. Shen, S. Long, M. Wu, H. Shih, Q. Zheng, N. Yen, C.
%       Tung, and H. Liu. The empirical mode decomposition and the Hilbert
%       spectrum for nonlinear and non-stationary time series analysis.
%       Proc. R. Soc. Lond. A 1998 454: 903-995
%   [3] G. Rilling, P. Flandrin, and P. Gonc??alv`es. On empirical mode
%       decomposition and its algorithms. IEEE-EURASIP workshop on
%       nonlinear signal and image processing NSIP-03, Grado (I), 2003.

% Thomas Oberlin
% 12.2011
% thomas.oberlin@imag.fr
%
%NOTES
% Modified starting 20-Feb-2014 by Daniel J. GIbson from emdos.m as it
% appeared in http://www-ljk.imag.fr/membres/Thomas.Oberlin/EMDOS.tar.gz,
% referenced in [1]. There is still some code pertaining to the variable
% 'liss' that is probably dead code, but I haven't verified that, so I left
% it in.
%EXAMPLES
% IMFCA1 = memd_emd(CA1theta, 'preprocess', @do_nothing_pre, ...
%   'postprocess', @do_nothing_post, 'pre_params', '~~barf~~');


% Gets the parameter
[s, stop, alpha, maxmodes, t, liss, postprocess, preprocess, ...
    pre_params, localEMDfunc, localEMDparam] = init(varargin{:});

k = 1;
r=s;
imf = [];
preprocess_auxdata = [];

%main loop : requires at least 3 extrema to proceed
while ~ memd_stop_emd(r) && (k < maxmodes+1 || maxmodes == 0)
    old_r = r;
    tmp = r;
    % padding
    [indmin,indmax] = memd_extr(tmp);
    indmin=indmin+floor(liss(k));
    indmax=indmax+floor(liss(k));
    
    % Sifting
    stop_sift=0;
    aux=0;
    
    if ~isempty(preprocess)
        [r, preprocess_auxdata] = feval(preprocess, r, pre_params);
    end
    
    while ~stop_sift
        [tmin,tmax,mmin,mmax] = memd_boundary_conditions(indmin,indmax,t,r,r,6);
        envmin = interp1(tmin,mmin,t,'spline');
        envmax = interp1(tmax,mmax,t,'spline');
        envmoy = (envmin+envmax)/2;
        if ~isempty(localEMDfunc)
            w = feval(localEMDfunc, r, envmax, envmin, alpha, ...
                localEMDparam);
            nr = r - w .* envmoy;
            stop_sift = all(w == 0);
        else
            nr = r - envmoy;
            switch(stop)
                case 'f'
                    % Flandrin
                    amp = mean(abs(envmax-envmin))/2;
                    sx = abs(envmoy)./amp;
                    stop_sift = ~(mean(sx > alpha) > 0.05 | any(sx > 10*alpha));
                case 'h'
                    % Huang
                    stop_sift = norm(nr-r)/(norm(r)+eps) < alpha;
            end
        end
        
        if ~stop_sift
            r=nr; % Replaces signal with new signal.
            aux=aux+1;
            [indmin,indmax] = memd_extr(r); % Replaces min and max indices with new ones.
        end
    end
    
    % Defining IMF, possibly with some postprocessing.
    if ~isempty(postprocess)
        imf(k,:) =  feval(postprocess, r, preprocess_auxdata); %#ok<AGROW>
    else
        imf(k,:) =  r; %#ok<AGROW>
    end
    r = old_r - r; % Defining signal as signal minus IMF.
    k = k+1;
    
end
ort = memd_io(s,imf);

% Residue
if ~isempty(postprocess)
    imf(k,:) =  feval(postprocess, r, preprocess_auxdata); 
else
    imf(k,:) =  r; 
end

end

function [s, stop, alpha, maxmodes, t, liss, postprocess, preprocess, ...
    pre_params, localEMDfunc, localEMDparam] = init(varargin)
    
% INIT : internal function for the initialization of the parameters.
% Returns <s> as a row vector.

s = varargin{1};
if nargin == 2
  if isstruct(varargin{2})
    inopts = varargin{2};
  else
    error('when using 2 arguments the first one is the analyzed signal and the second one is a struct object describing the options')
  end
elseif nargin > 2
  try
    inopts = struct(varargin{2:end});
  catch
    error('bad argument syntax')
  end
end

% Default parameters.

defopts.stop = 'f';
defopts.alpha = 0.05;
defopts.maxmodes = 20;
defopts.t = 1:max(size(s));
defopts.liss = 0;
defopts.postprocess = [];
defopts.preprocess = [];
defopts.pre_params = [];
defopts.localEMDfunc = [];
defopts.localEMDparam = {};
opt_fields = {'stop','alpha','maxmodes','t','liss','postprocess',...
    'preprocess','pre_params','localEMDfunc','localEMDparam'};
opts = defopts;

if(nargin==1)
  inopts = defopts;
elseif nargin == 0
  error('not enough arguments')
end

% Some checking
names = fieldnames(inopts);
for nom = names'
  if ~any(strcmpi(char(nom), opt_fields))
    error(['bad option field name: ',char(nom)]);
  end
  % Et modification des param??tres rentr??s
  if ~isempty(eval(['inopts.',char(nom)])) % empty values are discarded
    eval(['opts.', char(nom), ' = inopts.', char(nom),';']);
  end
end

% Mise ?? jour
stop = opts.stop;
alpha = opts.alpha;
maxmodes = opts.maxmodes;
t = opts.t;
liss = opts.liss;
postprocess = opts.postprocess;
preprocess = opts.preprocess;
pre_params = opts.pre_params;
localEMDfunc = opts.localEMDfunc;
localEMDparam = opts.localEMDparam;


%% Syntax check
% s
if ~isvector(s)
  error('The signal S must have only one row or one column')
end
if size(s,1) > 1
  s = s.';
end

% t
if ~isvector(t)
  error('option field T must have only one row or one column')
end
if ~isreal(t)
  error('time instants T must be a real vector')
end
if size(t,1) > 1
  t = t';
end
if (length(t)~=length(s))
  error('X and option field T must have the same length')
end

% liss
if max(size(liss))>1
    if size(liss,1)>1
        liss=liss';
    end
    if size(liss,2)<maxmodes
        liss = [liss liss(end)*ones(1,maxmodes - length(liss))];
    end
else
    liss = [liss(1) zeros(1,maxmodes - 1)];
end


end
