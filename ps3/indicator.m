function is = indicator(point, min, max, type_min, type_max)
% INDICATOR(point, min, max, type_min, type_max) is 1 if point is between
% min and max values. You have to control endpoints with 'closed' or 'open'
% arguments. Otherwise it raises an error.

is = 0;
if strcmp(type_min, 'closed')==1
    if strcmp(type_max, 'closed')==1
        if point >= min && point <= max
            is = 1;
        end
    elseif strcmp(type_max, 'open')==1
        if point >= min && point < max
            is = 1;
        end
    else
        error('Upper bound is not clear')
    end
    
elseif strcmp(type_min, 'open')==1
    if strcmp(type_max, 'closed')==1
       if point > min && point <= max 
           is = 1;
       end
    elseif strcmp(type_max, 'open')==1
        if point > min && point < max
            is = 1;
        end
    else
        error('Upper bound is not clear.')
    end
else
    error('Lower bound is not clear')
end

end
    