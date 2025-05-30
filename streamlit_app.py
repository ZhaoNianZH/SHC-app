
"""
STATE POINT ANALYSIS
by Tommy Thompson

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

Future Work
Ultimately, the goal would be to integrate operational data into this tool as the inputs, such that the
output graph and Q_RAS values update during set time intervals to inform operators of what the most
efficient pump rate should be.


"""

# Required Libraries #
import matplotlib.pyplot as plt
import numpy as np
# from mpmath import lambertw
from scipy.special import lambertw

import streamlit as st



##############


class StatepointAnalysis:
    """
    StatepointAnalysis 类用于执行二级处理的状态点分析。
    它通过计算污泥通量、溢流率和污泥回流率等参数，帮助操作员确定最高效的污泥回流泵速率。
    """

    def __init__(self, number_of_tanks, tank_area, test_type, sludge_volume_index, mixed_liquor_ss, q_in):
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
            "SVIGS": (0.245, 0.00296, 0.01073, 24.3)
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

    def plot_spa(self, q_in, q_ras, tank_area, number_of_tanks, sludge_volume_index, mixed_liquor_ss, solids_max=15):
        """
        绘制状态点分析图。

        参数:
            q_in (float): 进水流量（m³/h）。
            q_ras (float): 泥水回流流量（m³/h）。
            tank_area (float): 每个澄清器的面积（平方米）。
            number_of_tanks (int): 二级澄清器的数量。
            sludge_volume_index (float): 泥水体积指数（SVI）。
            mixed_liquor_ss (float): 混合液悬浮固体浓度（g/L）。
            solids_max (float): 最大污泥浓度（g/L），默认为 20。
        """
        x_overflow_rate, y_overflow_rate = self.surface_overflow_calc(q_in, tank_area, number_of_tanks, sludge_volume_index)
        x_underflow_rate, y_underflow_rate = self.solids_underflow_calc(q_in, q_ras, mixed_liquor_ss, tank_area, number_of_tanks)
        x_flux, y_flux = self.settle_flux(sludge_volume_index, solids_max)
        surface_overflow_rate = y_overflow_rate[1] / x_overflow_rate[1]
        state_point = surface_overflow_rate * mixed_liquor_ss

        # 创建一个新的 figure 对象
        fig, ax = plt.subplots()
        ax.plot(x_overflow_rate, y_overflow_rate, "g", label="Surface Overflow Rate")
        ax.plot(x_underflow_rate, y_underflow_rate, "orange", label="Solids Underflow Rate")
        ax.plot(x_flux, y_flux, "b", label="Settle Flux")
        ax.plot(mixed_liquor_ss, state_point, "ro", label="State Point")
        ax.axis([0, solids_max, 0, y_overflow_rate[1]])
        ax.grid(True)

        ras_sub = r'$Q_{RAS}$'
        txt = "%s = %s m3/h" % (ras_sub, q_ras)

        ax.set_title(f"State Point Analysis: {txt}")
        ax.set_xlabel("Solids Concentration (g/L)")
        ax.set_ylabel("Solids Flux (kg/m2/h)")
        #ax.text(0.5 * solids_max, 0.5 * y_overflow_rate[1], txt, size=12, backgroundcolor="white", bbox=dict(facecolor="white"))
        ax.legend()

        # 使用 st.pyplot 显示图表
        st.pyplot(fig)

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
        q_ras = float((self.Qin * self.MLSS * k_sp) / lamb_w.real)  # m3/h; convert to Python float
        op_qras = self.oper_round(self.Qin, q_ras)  # m3/h; convert to Python float

        # Plot State Point Analysis Graph
        underflow_x_intercept = self.solids_underflow_calc(self.Qin, q_ras, self.MLSS, self.TANK_AREA, self.NUMBER_OF_TANKS)[0][1] # g/L

        #solids_max = round(underflow_x_intercept / 5) * 5 + 2 
        #self.plot_spa(self.Qin, op_qras, self.TANK_AREA, self.NUMBER_OF_TANKS, self.SVI, self.MLSS, solids_max=solids_max) # g/L; Max solids conc. on x-axis rounded to multiple of 5
        #print(q_ras, op_qras)
        return q_ras, op_qras


# Streamlit 应用
def main():
    st.title("State Point Analysis Calculator")
    
    # 输入框
    number_of_tanks = st.number_input("Number of Secondary Clarifiers", value=3, min_value=1)
    tank_area = st.number_input("Area of Each Clarifier (m²)", value=1986.0, min_value=1.0)
    test_type = st.selectbox("Test Type", ["SVISN", "SVISS", "SVIGN", "SVIGS"])
    sludge_volume_index = st.number_input("Sludge Volume Index (SVI)", value=140.0, min_value=0.0)
    mixed_liquor_ss = st.number_input("Mixed Liquor Suspended Solids (g/L)", value=4.3, min_value=0.0)
    q_in = st.number_input("Influent Flow Rate (m³/h)", value=5678.0, min_value=0.0)

    # 计算按钮
    if st.button("Calculate"):
        # 创建 StatepointAnalysis 实例
        spa = StatepointAnalysis(number_of_tanks, tank_area, test_type, sludge_volume_index, mixed_liquor_ss, q_in)
        q_ras, op_qras = spa.find_state_point()
        
        # 显示结果
        st.write(f"Optimal Q_RAS: {op_qras} m³/h")
        
        # 绘制图形
        fig, ax = plt.subplots()
        spa.plot_spa(q_in, op_qras, tank_area, number_of_tanks, sludge_volume_index, mixed_liquor_ss)
        #st.pyplot(fig)

if __name__ == "__main__":
    main()
