function varargout = timing(varargin)
% TIMING  display begin time, end time, and calculate duration of a command
%

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 10/13/2005
% $Id: timing.m,v 1.2 2012/03/11 01:16:14 jonp Exp $
%**************************************************************************%
  
  t0 = tic();
  %t0 = clock;
    
  if ( nargin == 0 ),
    fprintf(1, '\nstart time: %s\n\n', datestr(now));

    if ( nargout >= 1 ),
      varargout{1} = t0;
    end;
    
    return;
  end;
    
  t0 = varargin{1};
  
  t1 = toc(t0);
  %t1 = clock;

  fprintf(1, '\n  end time: %s\n', datestr(now));
  
  runtime_seconds = t1;
  %runtime_seconds = etime(t1, t0);
  
  runtime_str = sprintf('[[ %02dh %02dm %02ds ]]', fix(runtime_seconds/60/60), ...
                        rem(fix(runtime_seconds/60), 60), ...
                        rem(fix(runtime_seconds), 60));


  if ( ~isempty(inputname(1)) ),
    fprintf(1, '\n\n  run time:  %s = %s\n\n', inputname(1), runtime_str);
  else,
    fprintf(1, '\n\n  run time:  %s\n\n', runtime_str);
  end;

  % flush diary
  if ( strcmp(get(0, 'Diary'), 'on') ),
    diary off; diary on;
  end;

  if ( nargout >= 1 ),
    varargout{1} = tic();
  end;


  return;


  %************************************************************************%
  %%% $Source: /usr/fscvsroot/dev/gradient_nonlin_unwarp/timing.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
