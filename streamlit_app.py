
"""
STATE POINT ANALYSIS

Background
Wastewater treatment is a biochemical and mechanical process in which municipal and commercial wastewater
undergoes multiple, consecutive stages of separation to remove suspended and dissolved solids before
discharging treated effluent into a natural waterway.  The process train typically consists of:
     - Racks and screens to capture trash and other floatables
     - Primary sedimentation to settle out heavier suspended sediment or waste solids (sludge)
     - Aeration to promote bacterial growth for consuming remaining solids
     - Secondary sedimentation (clarification) to settle out waste-bacteria aggregates
     - Discharge of treated effluent into a local waterbody for finishing by natural microorganisms

In most setups, the aeration-secondary treatment stage occur simultaneously in a continuous-flow loop
process.  That is, part of the sludge blanket that settles out at the bottom of the secondary clarifier
is intentionally pumped back in to the aeration tank to more efficiently promote the growth of the
necessary microorganisms.  This is known as "Return Activated Sludge (RAS)," referring to it being pumped
back into the previous stage (return) and primed full of beneficial microorganisms (activated).

The volumetric flow rate at which the RAS is pumped into the aeration tank is controlled by the
Wastewater Operator.  The ultimate goal is to maintain a stable yet ever-present sludge blanket at the
bottom of the secondary clarifier: pump too much and it disappears; pump too little and it accumulates
and eventually leaves the clarifier as effluent.

To more easily and efficiently maintain this preferred operating point, an idealized model for secondary
treatment was created known as "State Point Analysis."  Given parameters such as the solids concentration
of the aeration tank (Mixed Liquor Suspended Solids) and the sludge volume index (a measure of how readily
the solids settle in water), the reactor hydraulics of the clarifier can be plotted to visualize whether
the operator is running the RAS pump too low or too high.

Description
The following module aims to maximize energy efficiency, minimize risk of losing solids, and simplify the
operation of secondary treatment by analytically determining the precise point below which the sludge
blanket would begin to rise and threaten to exit the system.

STATE POINT ANALYSIS
Optimal RAS calculation by Tommy Thompson
Sludge Handling Classification by ZHZN
Streamlit implementation by ZHZN

Last edit 2025-06-08
"""

# Required Libraries #
import matplotlib.pyplot as plt
import numpy as np
from scipy.special import lambertw
import streamlit as st



##############


class StatepointAnalysis:
    """
    StatepointAnalysis 类用于执行二级处理的状态点分析。
    它通过计算污泥通量、溢流率和污泥回流率等参数，帮助操作员确定最高效的污泥回流泵速率。
    """

    def __init__(self, number_of_tanks:int, tank_area:float, test_type:str, sludge_volume_index:int, mixed_liquor_ss:float, q_in:float, q_ras:float = None):
        """
        初始化 StatepointAnalysis 类实例。

        参数:
            number_of_tanks (int): 二级澄清器的数量。
            tank_area (float): 每个澄清器的面积（平方米）。
            test_type (str): 测试类型，用于选择常数集。
            sludge_volume_index (float): 泥水体积指数（SVI）。
            mixed_liquor_ss (float): 混合液悬浮固体浓度（g/L）。
            q_in (float): 进水流量（m³/h）。
        """
        self.NUMBER_OF_TANKS = number_of_tanks  # 二级澄清器的数量
        self.TANK_AREA = tank_area  # 每个澄清器的面积（平方米）
        self.test_type = test_type  # 测试类型
        self.SVI = sludge_volume_index  # 泥水体积指数（SVI）
        self.MLSS = mixed_liquor_ss  # 混合液悬浮固体浓度（g/L）
        self.Qin = q_in  # 进水流量（m³/h）
        self.Qras = q_ras # 污泥回流流量（m³/h）
        self.alpha, self.beta, self.delta, self.gamma = self.def_constants(test_type)  # 常数集

    def def_constants(self, test_type="SVISN"):
        """
        根据测试类型定义常数集。

        参数:
            test_type (str): 测试类型，默认为 "SVISN"。

        返回:
            tuple: 包含四个常数的元组 (alpha, beta, delta, gamma)。
        """
        test_dict = {
            "SVISN": (0.261, 0.00170, 0.00370, 14.9),
            "SVISS": (0.211, 0.00236, 0.00593, 14.6),
            "SVIGN": (0.351, 0.00058, 0.00602, 18.2),
            "SVIGS": (0.245, 0.00296, 0.01073, 24.3),
            "Custom": (0.165, 0.00159, -1.871/self.SVI, 1)
        }
        return test_dict[test_type]

    def surface_overflow_calc(self, q_in, tank_area, number_of_tanks, sludge_volume_index):
        """
        计算表面溢流率。

        参数:
            q_in (float): 进水流量（m³/h）。
            tank_area (float): 每个澄清器的面积（平方米）。
            number_of_tanks (int): 二级澄清器的数量。
            sludge_volume_index (float): 泥水体积指数（SVI）。

        返回:
            tuple: 包含 x 和 y 坐标的元组 (x, y)。
        """
        
        #       1.25*γ*e^(-(δ*SVI + 1))
        # Y_2 = -----------------------
        #             α + β*SVI
        # 对Flux 求导，令其为0，即得到 solids_peak: x = 1 / (aα + β*SVI), 将 x 带入 Flux peak flux,为保证作图后 sor y轴上的
        # 范围是高于peak flux 的，对 peak flux * 1.25 得到 y_sor_2； 反代入sor 求得对应的 x_2;
        x_sor_1, y_sor_1 = 0, 0
        y_sor_2 = 1.25 * self.gamma * np.exp(-(self.delta * sludge_volume_index + 1)) / (self.alpha + self.beta * sludge_volume_index)
            # unit: kg/m2/h

        #       Y_2*N*A
        # X_2 = -------
        #        Q_in
        x_sor_2 = y_sor_2 * number_of_tanks * tank_area / q_in
        x = (x_sor_1, x_sor_2)
        y = (y_sor_1, y_sor_2)
        return x, y

    def solids_underflow_calc(self, q_in, q_ras, mixed_liquor_ss, tank_area, number_of_tanks):
        """
        计算污泥回流率。

        参数:
            q_in (float): 进水流量（m³/h）。
            q_ras (float): 泥水回流流量（m³/h）。
            mixed_liquor_ss (float): 混合液悬浮固体浓度（g/L）。
            tank_area (float): 每个澄清器的面积（平方米）。
            number_of_tanks (int): 二级澄清器的数量。

        返回:
            tuple: 包含 x 和 y 坐标的元组 (x, y)。
            x[1]: x-intercept of the line. 二沉池回流污泥浓度Xr。
            y[0]: y-intercept of the line. 二沉池表面固体负荷。
        """
        # Assign for x- and y-intercepts
        x_sur_1, y_sur_2 = 0, 0

        #       (Q_in + Q_ras)*MLSS
        # Y_1 = -------------------  SLR
        #               N*A
        y_sur_1 = (q_in + q_ras) * mixed_liquor_ss / (number_of_tanks * tank_area)
        # unit: kg/m2/h

        #       Y_1*N*A
        # X_2 = -------  X_return
        #        Q_ras
        x_sur_2 = y_sur_1 * number_of_tanks * tank_area / q_ras
        # unit: kg/m3

        x = (x_sur_1, x_sur_2)
        y = (y_sur_1, y_sur_2)

        return x, y
    def flux(self, solids_conc, sludge_volume_index):
        """
        计算污泥通量。

        参数:
            solids_conc (float): 污泥浓度（g/L）。
            sludge_volume_index (float): 泥水体积指数（SVI）。

        返回:
            float: 污泥通量（kg/m²/h）。
        """
        y = self.gamma * solids_conc * np.exp(-(self.delta * sludge_volume_index
                                               + (self.alpha + self.beta * sludge_volume_index) * solids_conc))
        return y

    def settle_flux(self, sludge_volume_index, solids_max=15):
        """
        计算污泥通量曲线。

        参数:
            sludge_volume_index (float): 泥水体积指数（SVI）。
            solids_max (float): 最大污泥浓度（g/L），默认为 20。

        返回:
            tuple: 包含 x 和 y 坐标的元组 (x, y)。
        """
        x = np.linspace(0, solids_max, num=100)
        y = self.flux(x, sludge_volume_index)
        return x, y

    def oper_round(self, ref, val):
        """
        对泥水回流流量进行四舍五入操作。

        参数:
            ref (float): 参考值，用于确定精度。
            val (float): 需要四舍五入的值。

        返回:
            float: 四舍五入后的值。
        """
        if ref < 10:
            precision = 0.10
        elif ref < 100:
            precision = 0.25
        else:
            precision = 0.50
        op_qras = np.ceil(val / precision) * precision
        return op_qras

    def find_state_point(self):
        """
        寻找状态点并绘制分析图。

        返回:
            tuple: 包含泥水回流流量 (q_ras) 和四舍五入后的泥水回流流量 (op_qras) 的元组。
        """
        # Solids concentration at which flux peaks (dy/dx = 0)
        flux_peak = 1 / (self.alpha + self.beta * self.SVI)

        # Two points to the right of peak for exponential approximation
        x_sp_1, x_sp_2 = map(lambda x: x * flux_peak, (3, 6)) # g/L; solids concentrations
        y_sp_1, y_sp_2 = map(self.flux, [x_sp_1, x_sp_2], [self.SVI] * 2) # kg/m2/h; corresponding fluxes

        # Y = A*e^(k*X) --> ln(Y) = k*X + lnA (linearization of exponential function)
        x_arr = np.array([[x_sp_1, 1], 
                          [x_sp_2, 1]])# k coefficients plus lnA constants
        y_arr = np.array([np.log(y_sp_1), 
                          np.log(y_sp_2)]) # Natural logarithms of Y values
        

        k_sp, nat_log_a = np.linalg.solve(x_arr, y_arr)  # L/g, ln(kg/m2/h); estimates for line
        a = np.exp(nat_log_a) # lkg/m2/h; exponential coefficient

        # Apply Lambert W function to calculate q_ras
        w_input = -self.Qin * self.MLSS / (self.NUMBER_OF_TANKS * self.TANK_AREA * a * np.exp(self.MLSS * k_sp + 1))
        lamb_w = min(map(lambertw, [w_input] * 2, [0, -1])) # W func gives two solutions, take the smaller

        # Calculate operating variable, q_ras
        q_ras_opt = float((self.Qin * self.MLSS * k_sp) / lamb_w.real)  # m3/h; convert to Python float
        op_qras = self.oper_round(self.Qin, q_ras_opt)  # m3/h; convert to Python float

        return q_ras_opt, op_qras

    def settle_flux_total(self, sludge_volume_index, solids_max=15):
        """
        二沉池的总固体通量，重力沉降和底流出流之和。

        参数:
            sludge_volume_index (float): 泥水体积指数（SVI）。
            solids_max (float): 最大污泥浓度（g/L），默认为 15。
            q_ras (float): 污泥回流流量 (m3/h)。
 
        返回:
            tuple: 包含 x 和 y 坐标的元组 (x, y)。
        """

        x = np.linspace(0, solids_max, num=100)
        y = self.flux(x, sludge_volume_index) + self.Qin / (self.TANK_AREA * self.NUMBER_OF_TANKS) * self.Qin / self.Qras * x 
        return x, y



    def plot_spa(self, q_in, q_ras, tank_area, number_of_tanks, sludge_volume_index, mixed_liquor_ss, solids_max=15):
        """
        绘制状态点分析图。
        """
        # MLSS Qin q_ras SVI SOR 
        # 计算SOR曲线
        x_overflow_rate, y_overflow_rate = self.surface_overflow_calc(q_in, tank_area, number_of_tanks, sludge_volume_index)
        solids_max = x_overflow_rate[1]*1.5

        # 计算UFO曲线
        x_underflow_rate, y_underflow_rate = self.solids_underflow_calc(q_in, q_ras, mixed_liquor_ss, tank_area, number_of_tanks)
        x_r = x_underflow_rate[1]  

        # 计算重力沉降曲线 
        x_flux, y_flux = self.settle_flux(sludge_volume_index, solids_max)
        
        # 计算状态点
        surface_overflow_rate = y_overflow_rate[1] / x_overflow_rate[1]
        state_point = surface_overflow_rate * mixed_liquor_ss

        # 创建一个新的 figure 对象
        fig, ax = plt.subplots()
        ax.plot(x_overflow_rate, y_overflow_rate, "g", label="Surface Overflow Rate")
        ax.plot(x_underflow_rate, y_underflow_rate, "orange", label="Solids Underflow Rate")
        ax.plot(x_flux, y_flux, "b", label="Settle Flux")
        ax.plot(mixed_liquor_ss, state_point, "ro", label="State Point", markersize=5)
        ax.vlines(mixed_liquor_ss, 0, y_underflow_rate[0]*1.1, "r","--")
        ax.axis([0, solids_max, 0, y_underflow_rate[0]*1.1])
        ax.grid(True)

        # 固定图例位置到右上角
        ax.legend(loc='upper right', fontsize=8)  # 修改：缩小图例字体大小

        # 调整 X 轴和 Y 轴刻度为 1
        ax.xaxis.set_major_locator(plt.MultipleLocator(2))
        ax.yaxis.set_major_locator(plt.MultipleLocator(1))

        txt = f"SLR: {y_underflow_rate[0]:.1f} (kg/m$^2$/h), SOR: {surface_overflow_rate:.2f} m/h, RAS: {q_ras/q_in:.1f}, Q$_r$: {q_ras:.1f} m$^3$/h, X$_r$: {x_r:.1f}g/L"

        ax.set_title(f"{txt}", fontsize=10)  # 修改：缩小标题字体大小
        ax.set_xlabel("Solids Concentration (g/L)")
        ax.set_ylabel("Solids Flux (kg/m$^2$/h)")

        return fig,txt



# Streamlit 应用
def page_ras_optimal_calculation():
    """
    RAS Optimal Calculation 页面。
    计算RAS最优解。
  
    """
    st.title("Optimal RAS Calculation")
    
    # 输入框，增大步长
    number_of_tanks = st.number_input("Number of Secondary Clarifiers", value=3, min_value=1, step=1)
    tank_area = st.number_input("Area of Each Clarifier (m²)", value=1986.0, min_value=1.0, step=100.0)
    test_type = st.selectbox("Test Type", ["SVISN", "SVISS", "SVIGN", "SVIGS","Custom"])
    sludge_volume_index = st.number_input("Sludge Volume Index (SVI)", value=140.0, min_value=0.0, step=5.0)
    mixed_liquor_ss = st.number_input("Mixed Liquor Suspended Solids (g/L)", value=4.3, min_value=0.0, step=0.1)
    q_in = st.number_input("Influent Flow Rate (m³/h)", value=5678.0, min_value=0.0, step=100.0)

    # 使用 session_state 存储输入参数和结果
    if "ras_params" not in st.session_state:
        st.session_state.ras_params = {
            "number_of_tanks": number_of_tanks,
            "tank_area": tank_area,
            "test_type": test_type,
            "sludge_volume_index": sludge_volume_index,
            "mixed_liquor_ss": mixed_liquor_ss,
            "q_in": q_in,
            "result": None
        }

    # 检查参数是否发生变化
    params_changed = any([
        number_of_tanks != st.session_state.ras_params["number_of_tanks"],
        tank_area != st.session_state.ras_params["tank_area"],
        test_type != st.session_state.ras_params["test_type"],
        sludge_volume_index != st.session_state.ras_params["sludge_volume_index"],
        mixed_liquor_ss != st.session_state.ras_params["mixed_liquor_ss"],
        q_in != st.session_state.ras_params["q_in"]
    ])

    # 计算按钮
    if st.button("Calculate", key=f"calculate_button_{number_of_tanks}_{tank_area}_{test_type}_{sludge_volume_index}_{mixed_liquor_ss}_{q_in}") or params_changed:
        # 更新 session_state 中的参数
        st.session_state.ras_params = {
            "number_of_tanks": number_of_tanks,
            "tank_area": tank_area,
            "test_type": test_type,
            "sludge_volume_index": sludge_volume_index,
            "mixed_liquor_ss": mixed_liquor_ss,
            "q_in": q_in
        }

        # 创建 StatepointAnalysis 实例
        spa = StatepointAnalysis(
            number_of_tanks, tank_area, test_type, sludge_volume_index, mixed_liquor_ss, q_in
        )
        q_ras_opt, op_qras = spa.find_state_point()
        fig,txt = spa.plot_spa(q_in, op_qras, tank_area, number_of_tanks, sludge_volume_index, mixed_liquor_ss)
        # 更新结果
        st.session_state.ras_params["result"] = {
            "q_ras_opt": q_ras_opt,
            "op_qras": op_qras,
            "fig": fig
        }

    # 显示结果
    if st.session_state.ras_params["result"]:
        st.write(f"Optimal Q_RAS: {st.session_state.ras_params['result']['op_qras']} m³/h", font_size=10)  # 修改：缩小文本字体大小
        st.pyplot(st.session_state.ras_params["result"]["fig"])


def page_shc_analysis():
    """
    SHC Analysis 页面
    计算当前工况下的状态点，并根据SHC状态图进行评价。
    """
    st.title("SHC Analysis")
    
    # 输入框，增大步长
    number_of_tanks = st.number_input("Number of Secondary Clarifiers", value=3, min_value=1, step=1)
    tank_area = st.number_input("Area of Each Clarifier (m²)", value=1986.0, min_value=1.0, step=100.0)
    test_type = st.selectbox("Test Type", ["SVISN", "SVISS", "SVIGN", "SVIGS","Custom"])
    sludge_volume_index = st.number_input("Sludge Volume Index (SVI)", value=140.0, min_value=0.0, step=5.0)
    mixed_liquor_ss = st.number_input("Mixed Liquor Suspended Solids (g/L)", value=4.3, min_value=0.0, step=0.1)
    q_in = st.number_input("Influent Flow Rate (m³/h)", value=5678.0, min_value=0.0, step=100.0)
    q_ras = st.number_input("Underflow Rate (m³/h)", value=5678.0, min_value=0.0, step=100.0)

    # 使用 session_state 存储输入参数和结果
    if "shc_params" not in st.session_state:
        st.session_state.shc_params = {
            "number_of_tanks": number_of_tanks,
            "tank_area": tank_area,
            "test_type": test_type,
            "sludge_volume_index": sludge_volume_index,
            "mixed_liquor_ss": mixed_liquor_ss,
            "q_in": q_in,
            "q_ras": q_ras,
            "result": None
        }

    # 检查参数是否发生变化
    params_changed = any([
        number_of_tanks != st.session_state.shc_params["number_of_tanks"],
        tank_area != st.session_state.shc_params["tank_area"],
        test_type != st.session_state.shc_params["test_type"],
        sludge_volume_index != st.session_state.shc_params["sludge_volume_index"],
        mixed_liquor_ss != st.session_state.shc_params["mixed_liquor_ss"],
        q_in != st.session_state.shc_params["q_in"],
        q_ras != st.session_state.shc_params["q_ras"]
    ])

    # 计算按钮
    if st.button("Calculate", key=f"calculate_button_{number_of_tanks}_{tank_area}_{test_type}_{sludge_volume_index}_{mixed_liquor_ss}_{q_in}_{q_ras}") or params_changed:
        # 更新 session_state 中的参数
        st.session_state.shc_params = {
            "number_of_tanks": number_of_tanks,
            "tank_area": tank_area,
            "test_type": test_type,
            "sludge_volume_index": sludge_volume_index,
            "mixed_liquor_ss": mixed_liquor_ss,
            "q_in": q_in,
            "q_ras": q_ras
        }

        # 创建 StatepointAnalysis 实例
        spa = StatepointAnalysis(
            number_of_tanks, tank_area, test_type, sludge_volume_index, mixed_liquor_ss, q_in, q_ras
        )
        fig,txt = spa.plot_spa(q_in, q_ras, tank_area, number_of_tanks, sludge_volume_index, mixed_liquor_ss)

        # 更新结果
        st.session_state.shc_params["result"] = {
            "fig": fig
        }

    # 显示结果
    if st.session_state.shc_params["result"]:
        st.pyplot(st.session_state.shc_params["result"]["fig"])
        #st.write(txt)

    # 添加图片
    st.image("pictures/SHC_chart.png", caption="SHC Chart", use_container_width=True)
    st.write("**SHC I (underflow condition):**  if the underflow line happens to fall outside of the flux curve, we cannot get the solids out from the clarifier (settling is insufficient)")    
    st.write("**SHC II (feed load condition):**  if the state point happens to fall outside of the flux curve, the load in the feed layer is more than what the clarifier can handle at that point (so it goes up eventually, as the clarifier is overloaded).")


def main():
    # 通过左侧导航栏切换页面
    page = st.sidebar.radio("Select Page", [ "SHC analysis","Optimal RAS calculation"])
    
    if page == "Optimal RAS calculation":
        page_ras_optimal_calculation()
    elif page == "SHC analysis":
        page_shc_analysis()

if __name__ == "__main__":
    main()
