function hh = titlef(str)

if (isnumeric(str))
    str = sprintf('Iteration: %d',str);
end
hh = title(str,'FontSize',14);