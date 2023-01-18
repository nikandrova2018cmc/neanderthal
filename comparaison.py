

def comparaison(reality,calculated):
    total_gen = sum(list(map(sum,reality)))
    errors = []
    not_inside = 0
    not_inside1 = 0
    for real in calculated:
        flag = True
        for calc in reality:
            start_real,end_real,start_calc,end_calc = real[0],real[1],calc[0],calc[1]
            maximum_margin = int(end_real-start_real)*0.75
            if abs(start_real-start_calc) < maximum_margin and abs(end_real-end_calc) < maximum_margin:
                error = abs(start_real-start_calc) + abs(end_real-end_calc)
                errors.append(error)
                flag = False
                break
            error = abs(start_real-end_real)
            errors.append(error)
        if flag:
            not_inside += 1
        not_inside1  +=1
    print(not_inside/not_inside1)
    error_rate = (sum(errors)/total_gen)*100
    return error_rate

