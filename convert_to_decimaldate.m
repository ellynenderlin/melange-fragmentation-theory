function [decidate] = convert_to_decimaldate(varargin)
%%% Converts input datestring to a decimal date. The default input format is
%%% a datestring in YYYYMMDD format. Another datestring format can be
%%% specified as long as it uses the same characters as Matlab's defaults. 

%days of year 
modays_norm = [31 28 31 30 31 30 31 31 30 31 30 31];
cumdays_norm = cumsum(modays_norm); cumdays_norm = [0 cumdays_norm(1:11)];
modays_leap = [31 29 31 30 31 30 31 31 30 31 30 31];
cumdays_leap = cumsum(modays_leap); cumdays_leap = [0 cumdays_leap(1:11)];

%assign inout arguments to variable names
datestring = varargin{1};
if nargin == 2
    format = varargin{2};
end

%check that the input date is a string
if ~ischar(datestring)
    error('Input date must be a char'); %exit function if the format is wrong
end

%check if an alternate datestring format is provided
if nargin == 2
    if isempty(strmatch(format,'YYYYMMDD'))
        %find starting index for the year
        year_start = strfind(format,'YYYY');
        if isempty(year_start)
            error('incorrect date input format: datestring must include YYYY for year'); 
        end
        
        %find starting index for the month
        mo_start = strfind(format,'MM');
        if isempty(mo_start)
            error('incorrect date input format: datestring must include MM for month'); 
        end
        
        %find starting index for the day
        day_start = strfind(format,'DD');
        if isempty(day_start)
            error('incorrect date input format: datestring must include DD for day'); 
        end
        
        %rearrange to match the desired input format
        newdate = [datestring(year_start:year_start+3),datestring(mo_start:mo_start+1),datestring(day_start:day_start+1)];
        clear datestring; datestring = newdate; format = 'YYYYMMDD'; clear newdate;
    end
end

%convert to decimal date
if mod(str2num(datestring(1:4)),4) == 0
    decidate = str2num(datestring(1:4))+((cumdays_leap(str2num(datestring(5:6)))+str2num(datestring(7:8)))/sum(modays_leap));
else
    decidate = str2num(datestring(1:4))+((cumdays_norm(str2num(datestring(5:6)))+str2num(datestring(7:8)))/sum(modays_norm));
end

end