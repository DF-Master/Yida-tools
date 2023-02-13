import streamlit as st
import numpy as np
from scipy import interpolate
import time


def interpolate_list_input(list_input,
                           x_list_pos=0,
                           interpolate_count=10,
                           mode="slinear",
                           int_out=True):
    #"nearest","zero"为阶梯插值
    #slinear 线性插值
    #"quadratic","cubic" 为2阶、3阶B样条曲线插值
    # list form:[x_list,y1_list,y2_list,……]

    x_list_np = list_to_numpy(list_input[x_list_pos])
    count_list_np = np.linspace(1, len(x_list_np), len(x_list_np))
    count_list_np_new = np.linspace(
        1, len(x_list_np),
        interpolate_count * len(x_list_np) - interpolate_count + 1)
    # Expend X list
    f_x = interpolate.interp1d(count_list_np, x_list_np, kind="slinear")
    x_list_new_np = f_x(count_list_np_new)

    del list_input[x_list_pos]
    y_all_list = []
    for i in range(len(list_input)):
        f = interpolate.interp1d(x_list_np,
                                 list_to_numpy(list_input[i]),
                                 kind=mode)
        y_all_list.append(list(f(x_list_new_np)))
    if x_list_pos < len(y_all_list):
        y_all_list.insert(x_list_pos, list(x_list_new_np))
    else:
        y_all_list.append(list(x_list_new_np))

    if int_out == True:
        # Round not int
        y_all_list = [[str(round(j)) for j in i] for i in y_all_list]
        return y_all_list

    else:
        return y_all_list


def Transform_list(l):
    return list(map(list, zip(*l)))


def list_to_numpy(list_input):
    return np.array(list_input)


if __name__ == '__main__':
    st.title("Converter-WebUI")
    st.caption('made by [Yida](https://github.com/DF-Master) --221101 update!',
               unsafe_allow_html=True)
    output = st.empty()

    st.subheader("List Transformer")
    transform_list_input = st.text_input(
        'Transform_List_Input(a,b,c,……;A,B,C,……)',
        value='1,11,21;2,22,42;3,33,63')
    if st.button('Transform'):
        try:

            st.write("###########################################")
            st.write(time.strftime("%a %b %d %H:%M:%S %Y", time.localtime()))
            st.write("###########################################")

            input_list = [
                i.split(",") for i in transform_list_input.split(";")
            ]

            output_list = Transform_list(input_list)

            output_str = ";".join([",".join(i) for i in output_list])

            st.write(output_str)
            output.markdown(' Transform Finished ')
        except:
            output.markdown(' Transform failed ')

    st.subheader("Interpolate Converter")

    base_set = st.text_input(
        'Interpolater(a,b,c,……;A,B,C,……)[A:x axis;B:y1 axis;C y2 axis]',
        value='1,11,21;2,22,42;3,33,63')
    x_axis_pos_set = st.text_input("X_axis_Pos,Default=2(Third)", value=2)
    interpolate_count_set = st.text_input("Interpolate_Count_Set", value=10)
    mode_set = st.radio("Mode_set",
                        ("slinear", "nearest", "zero", "quadratic", "cubic"))
    int_out_set = st.radio("Int_out_set", ("1", "0"))

    if st.button('Interpolate'):
        try:

            st.write("###########################################")
            st.write(time.strftime("%a %b %d %H:%M:%S %Y", time.localtime()))
            st.write("###########################################")

            input_list = [i.split(",") for i in base_set.split(";")]
            input_list = [[float(j) for j in i] for i in input_list]
            output_list = interpolate_list_input(
                input_list,
                x_list_pos=int(x_axis_pos_set),
                interpolate_count=int(interpolate_count_set),
                mode=mode_set,
                int_out=int(int_out_set))

            output_str = ";".join([",".join(i) for i in output_list])

            st.write(output_str)
            output.markdown(' Interpolate Finished ')
        except:
            output.markdown(' Interpolate failed ')

    # print(interpolate_list_input([[1, 11], [2, 22]]))

    st.subheader("List Difference[Keep First One]")

    diff_list_input = st.text_input('Diff_List_Input(a,b,c,……;A,B,C,……)',
                                    value='1,11,21;2,22,42;3,33,63')
    target_list_set = st.text_input('Target_List_Set(a,b,c,……)', value='1,3')
    if st.button('Diff'):
        try:

            st.write("###########################################")
            st.write(time.strftime("%a %b %d %H:%M:%S %Y", time.localtime()))
            st.write("###########################################")

            target_num_list = [int(i) for i in target_list_set.split(",")]
            input_list_all = diff_list_input.split(";")
            input_list = [[int(j) for j in input_list_all[i].split(",")]
                          for i in target_num_list]

            output_list = [[str(i[0])] +
                           [str(i[j + 1] - i[j]) for j in range(len(i) - 1)]
                           for i in input_list]

            for i in range(len(target_num_list)):
                input_list_all.insert(target_num_list[i],
                                      ",".join(output_list[i]))
                del input_list_all[target_num_list[i] + 1]
            output_str = ";".join(input_list_all)

            st.write(output_str)
            output.markdown(' List Diff Finished ')
        except:
            output.markdown(' List Diff failed ')

    st.subheader("List UV-3D-Control-Converter")

    simple_input = st.text_input('Simple_List_Input(a,b,c,……;A,B,C,……)',
                                 value='1,11,21;2,22,42;3,33,63')
    if st.button('UV-3D-Control-Convert'):
        try:

            st.write("###########################################")
            st.write(time.strftime("%a %b %d %H:%M:%S %Y", time.localtime()))
            st.write("###########################################")

            input_list = [i.split(",") for i in simple_input.split(";")]

            input_list = Transform_list(input_list)
            input_list = [[float(j) for j in i] for i in input_list]

            input_list = interpolate_list_input(
                input_list,
                x_list_pos=int(x_axis_pos_set),
                interpolate_count=int(interpolate_count_set),
                mode=mode_set,
                int_out=int(int_out_set))

            input_list = [[str(j) for j in i] for i in input_list]
            input_list = [",".join(i) for i in input_list]

            target_num_list = [int(i) for i in target_list_set.split(",")]
            input_list_all = input_list
            input_list = [[int(j) for j in input_list_all[i].split(",")]
                          for i in target_num_list]

            output_list = [[str(i[0])] +
                           [str(i[j + 1] - i[j]) for j in range(len(i) - 1)]
                           for i in input_list]

            for i in range(len(target_num_list)):
                input_list_all.insert(target_num_list[i],
                                      ",".join(output_list[i]))
                del input_list_all[target_num_list[i] + 1]
            output_list = Transform_list(
                [i.split(",") for i in input_list_all])
            print('here')
            output_str = ";".join([",".join(i) for i in output_list])
            target_num_list = [int(i) for i in target_list_set.split(",")]

            st.write(output_str)
            output.markdown(' UV-3D-Control-Convert Finished ')
        except:
            output.markdown(' UV-3D-Control-Convert failed ')
