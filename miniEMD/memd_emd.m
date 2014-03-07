function [imf,ort] = memd_emd(varargin)
% EMDOS : Computes the EMD by Optimisation on Splines, implements the
% method described in [1]. Uses either the standard EMD [2,3] or the OS
% algorithm [1].
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
%       maxmodes    <8>     : Maximum number of IMFs; <0> means no maximum.
%       postprocess <empty> : Handle to a function that accepts as
%                       arguments (1) the newly extracted IMF <r> and (2)
%                       the structure of auxiliary values returned by its
%                       matching preprocess function.  If there is no
%                       preprocessing function, the second argument will be
%                       empty. The function is called immediately following
%                       the termination of sifting for each IMF. <r> is a
%                       column vector when passed in to the postprocess
%                       function, which returns a modified version as a
%                       column vector. This is where you would e.g. remove
%                       the masking signal added by <preprocess>.
%       preprocess  <empty> : Handle to a function that accepts as
%                       arguments the residual signal <r> and a value
%                       (which could be a struct, scalar, cell array, or
%                       whatever) containing any necessary parameters. It
%                       must return a modified version of <r> that is used
%                       in place of the original <r>, and a value
%                       containing any auxiliary data that will be needed
%                       by the matching 'postprocess'.  The preprocess
%                       function is executed immediately before beginning
%                       the sifting for each IMF.  <r> is a column vector
%                       when passed to the preprocess function, and the
%                       modified version returned must also be a column
%                       vector. This is where you would e.g. add a masking
%                       signal.
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
% Modified 20-Feb-2014 by Daniel J. GIbson from emdos.m as it appeared in
% http://www-ljk.imag.fr/membres/Thomas.Oberlin/EMDOS.tar.gz, referenced in
% [1].
% There is still some code pertaining to the variable 'liss' that is
% probably dead code, but I haven't verified that, so I left it in.
%EXAMPLES
% IMFCA1 = memd_emd(CA1theta, 'preprocess', @do_nothing_pre, ...
%   'postprocess', @do_nothing_post, 'pre_params', '~~barf~~');


% Gets the parameter
[s,stop,alpha,maxmodes,t,liss,postprocess,preprocess,pre_params] = ...
    init(varargin{:});

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
        envmin = interp1(tmin,mmin,t,'spline'); % Creates min envelope.
        envmax = interp1(tmax,mmax,t,'spline'); % Creates max envelope.
        envmoy = (envmin+envmax)/2; % Creates mean of min and max envelopes (underlying nonstationary/possibly oscillatory trend).
        nr = r-envmoy; % nr stands for "new r".
        
        
        switch(stop) % Checking whether to stop sifting, using one of two criteria.
            case 'f'
                % Flandrin
                amp = mean(abs(envmax-envmin))/2; % Half of mean difference of max. and min. envelopes is the mean amp. of the signal.
                sx = abs(envmoy)./amp; % Divide underlying trend by this mean amplitude at each point.
                % Stop sifting if trend/amp is greater than alpha at fewer
                % than 5% of timepoints (i.e. trend/amp <= alpha 95% of the
                % time), and trend/amp is never bigger than 10*alpha. 
                stop_sift = ~(mean(sx > alpha) > 0.05 | any(sx > 10*alpha)); 
            case 'h'
                % Huang
                stop_sift = norm(nr-r)/(norm(r)+eps) < alpha; 
                % Stop sifting if the amplitude (2-norm) of the residual is
                % a fraction less than alpha of the signal started with.
        end
        
        if ~stop_sift
            r=nr; % Replaces signal with new signal.
            aux=aux+1;
            [indmin,indmax] = memd_extr(r); % Replaces min and max indices with new ones.
        end
    end
    
    % Defining IMF, possibly with some postprocessing.
    if ~isempty(postprocess)
        imf(k,:) =  feval(postprocess, r, preprocess_auxdata)'; %#ok<AGROW>
    else
        imf(k,:) =  r'; %#ok<AGROW>
    end
    r = old_r - r; % Defining signal as signal minus IMF.
    k = k+1;
    
end
ort = memd_io(s,imf);

% Residue
imf(k,:) = r';

end


function [s,stop,alpha,maxmodes,t,liss,postprocess,preprocess, ...
    pre_params] = init(varargin)
% INIT : internal function for the initialization of the parameters.


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
defopts.maxmodes = 8;
defopts.t = 1:max(size(s));
defopts.liss = 0;
defopts.postprocess = [];
defopts.preprocess = [];
defopts.pre_params = [];
opt_fields = {'stop','alpha','maxmodes','t','liss','postprocess',...
    'preprocess','pre_params'};
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
    error(['bad option field name: ',char(nom)])
  end
  % Et modification des param??tres rentr??s
  if ~isempty(eval(['inopts.',char(nom)])) % empty values are discarded
    eval(['opts.',lower(char(nom)),' = inopts.',char(nom),';'])
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
