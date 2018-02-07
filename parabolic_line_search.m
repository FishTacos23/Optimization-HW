function [a_star] = parabolic_line_search(a, f_a)
    [val, index2] = min(f_a);
    a_2 = a(index2);
    f_2 = f_a(index2);

    a_less = a;
    a_less(a_less<=a_2)=Inf;
    [val, index3] = min(a_less);
    a_3 = a(index3);
    f_3 = f_a(index3);

    a_more = a;
    a_more(a_more>=a_2)=0;
    [val, index1] = max(a_more);
    a_1 = a(index1);
    f_1 = f_a(index1);

    a_star = (f_1*(a_2^2-a_3^2)+f_2*(a_3^2-a_1^2)+f_3*(a_1^2-a_2^2))/(2*(f_1*(a_2-a_3)+f_2*(a_3-a_1)+f_3*(a_1-a_2)));
end