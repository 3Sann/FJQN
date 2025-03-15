function Kequal = k_equal(k,n)
    statusLeft = 1;
    statusRight = 1;
    for i=1:n
        statusLeft = statusLeft*nchoosek(2+k,2);
        statusRight = statusRight*nchoosek(2+k-1,2);
    end
    statusInFJ = statusLeft - statusRight;
    M = 2*n;
    statusAboutk = nchoosek(M - 1 + 1, M - 1);
    mindis = abs(statusInFJ - statusAboutk);
    Kequal = 1;
    for i = 2:1000
        statusAboutk = nchoosek(M - 1 + i, M - 1);
        dis = abs(statusInFJ - statusAboutk);
        if dis > mindis
            break
        end
        Kequal = i;
        if (Kequal - k)/k > 0.2       
            break
        end
        mindis = dis;
    end
end
