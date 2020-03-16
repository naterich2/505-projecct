#' ## Project
#' #### By Nathan Richman
#'Velocity profiles, Flow rate, and Wall shear stress can be calculated
#'for pulsatile flow by separating the steady and oscillatory profiles:
#'
#'Steady
#'======
#'
#'#### Velocity
#'
#'For steady flow in a tube we can make the following assumptions:
#'
#'1. $v_z = f(r)$
#'2. $v_r = v_{\theta} = 0$
#'3. Constant viscosity (newtonian), and density
#'4. $p = f(z)$
#'5. $g_z = 0$
#'
#'Under these assumptions, the $z$ component of the cylindrical
#'Navier-Stokes equations can be reduced and solved with the following
#'boundary conditions:
#'
#'1.  $\displaystyle\left.\frac{d v_z}{d r}\right|_{r=0} \ne \infty$
#'2.  $v_z(r = a) = 0$
#'$$\begin{aligned}
#'    0 &= -     \frac{d p}{dz} + \mu \left( \frac{1}{r}     \frac{d }{d r} \left(r     \frac{d v_z}{d r}\right)\right)\\
#'    \frac{r}{\mu} \frac{d p}{d z} &= \frac{d}{d r} \left( r \frac{d v_z}{d r} \right)\\
#'    \frac{d v_z}{d r} &=\frac{r}{2 \mu} \frac{d p}{d z} + \frac{c}{r} \rightarrow \textrm{BC 1: $c = 0$} \\
#'    v_z &= \frac{r^2}{4 \mu}     \frac{d p}{dz} + c, \ \textrm{BC 2: $v_z(r=a) = 0$}\\
#'    c  &= -\frac{a^2}{4 \mu}     \frac{d p}{dz} \\\end{aligned}$$
#'
#'#### Flow Rate
#'
#'Given the above velocity profile, we can calculate the flow rate $Q$
#'with: $$\begin{aligned}
#'    Q &= \int_0^a v_z(r) 2 \pi r dr\\
#'     &= \int_0^a \frac{a^2}{4 \mu} \frac{d p}{dz} \left(1 - \left(\frac{r}{a}\right)^2 \right) 2 \pi r dr\\
#'     &= \frac{2 \pi a^2}{4 \mu}\frac{d p}{dz}\int_0^a r - \frac{r^3}{a^2} dr\\
#'     &= \frac{2 \pi a^2}{4 \mu}\frac{d p}{dz}\left[ \frac{r^2}{2} - \frac{r^4}{4 a^2} \right]_0^a\\
#'     &=\frac{2 \pi a^2}{4 \mu}\frac{d p}{dz}\left( \frac{a^2}{2} - \frac{a^4}{4 a^2} \right)\\
#'     &=\frac{2 \pi a^4}{8 \mu}\frac{d p}{dz}\left(1 - 1/2 \right) \\
#'     &=\frac{ \pi a^4}{8 \mu}\frac{d p}{dz}\end{aligned}$$
#'
#'#### Shear Stress
#'
#'We can similarly calculate shear stress from the velocity profile:
#'$$\begin{aligned}
#'    \tau_{rz} &=2 \mu D_{rz} = 2 \mu * \frac{1}{2}\left(    \frac{\partial v_z}{\partial r}   +  \frac{\partial v_r}{dz}\right)\\
#'            &= \mu \left(\frac{\partial v_z(r,t)}{\partial r}\right)\\
#'              &= \frac{a^2}{4}\frac{d p}{dz} \frac{d}{d r} \left(1 - \frac{r^2}{a^2}\right) \\
#'              &= \frac{a^2}{4}     \frac{d p}{dz} \left( - \frac{2 r}{a^2}\right)\\
#'    \left.\tau_{rz}\right|_{r = a} &= \frac{a^2}{4}     \frac{d p}{dz} \left( - \frac{2 a}{a^2}\right)\\
#'                                   &= - \frac{a}{2}     \frac{d p}{dz}\end{aligned}$$
#'
#'Oscillatory
#'===========
#'
#'#### Velocity
#'
#'For unsteady flow in a tube, we can make the following assumptions:
#'
#'1.$v_z = f(r,t)$
#'2.$v_r = v_{\theta} = 0$
#'3.Constant viscosity (newtonian), and density
#'4.$p = f(z,t)$
#'5.$g_z = 0$
#'
#'These assumptions give us the following simplified $z$ component of the
#'cylindrical Navier-Stokes equations:
#'$$\rho \frac{\partial v_z}{\partial t}  = -     \frac{\partial p}{\partial z} + \mu \left[    \frac{1}{r}\frac{\partial}{\partial r}\left(r     \frac{\partial v_z}{\partial r}\right)\right]$$
#'If we assume the pressure gradient is a sum of sine or cosine waves and
#'represent this as a sum of exponentials:
#'$\frac{\partial p(z,t)}{\partial z} = p_s + p_o \sum_{n=1}^N (\cos n \omega t + i \sin n \omega t) = p_s + p_o \sum_{n=1}^N e^{i n \omega t}$,
#'where n is the harmonic of the wave. Then we can solve the governing
#'equation separately for each harmonic and sum the terms:
#'$$\begin{aligned}
#'    \rho \frac{\partial v_z}{\partial t}  &= -     p_o e^{i n \omega t} + \mu \left[    \frac{1}{r}\frac{\partial}{\partial r}\left(r     \frac{\partial v_z}{\partial r}\right)\right]\end{aligned}$$
#'Postulating that the velocity profile $v_z(r,t)$ will follow a similar
#'profile as the pressure gradient allows us to separate the $r$ and $t$
#'components: $v_z(r,t) = V(r) e^{i n \omega t}$. Using this separation of
#'variables we have that:
#'
#'-   $\frac{\partial v_z}{\partial t} = i n \omega V(r) e^{i n \omega t}$
#'-   $\frac{\partial v_z}{\partial r} = \frac{d V(r)}{d r} e^{i n \omega t}$
#'-   $\frac{d}{d r} r \frac{d V(r)}{r} = \frac{d V(r)}{d r} + r \frac{d^2 V(r)}{d r}$
#'
#'$$\begin{aligned}
#'    \rho i n \omega V(r) e^{i n \omega t}  &= -     p_o e^{i n \omega t} + \mu   \left[  \frac{d V(r)}{d r} + r \frac{d^2 V(r)}{d r}   \right] e^{i n \omega t} \\
#'    \frac{p_0}{\mu} &=\frac{d^2 V(r)}{d r^2} + \frac{1}{r} \frac{d V(r)}{d r} - \frac{\rho n i \omega V(r)}{\mu} \end{aligned}$$
#'For $\Omega_n = a \sqrt{\frac{\rho n \omega}{\mu}}$, this becomes:
#'$$\frac{d^2 V(r)}{d r^2} + \frac{1}{r} \frac{d V(r)}{d r} - \frac{i \Omega_n^2 V(r)}{a^2} = \frac{p_0}{\mu}$$
#'If we do a change of variables $\zeta = \Lambda \frac{r}{a}$, then
#'$\frac{d}{d r} = \frac{d}{d \zeta} \frac{d \zeta}{d r} = \frac{d}{d \zeta}\frac{\Lambda}{a}$,
#'and $\frac{d^2}{d r^2} = \frac{d^2}{d \zeta^2} \frac{\Lambda^2}{a^2}$.
#'If $\Lambda = \left(\frac{i - 1}{\sqrt{2}}\right) \Omega$, then we get:
#'$$\begin{aligned}
#'    \frac{\Lambda^2}{a^2} \frac{d^2 V(\zeta)}{d \zeta^2} + \frac{\Lambda}{a \zeta}\frac{\Lambda}{a} \frac{d V(\zeta)}{d \zeta} - \frac{i \Omega_n^2 V(\zeta)}{a^2} &= \frac{p_0}{\mu}\\
#'    \frac{-i \Omega_n^2}{a^2} \frac{d^2 V(\zeta)}{d \zeta^2} + \frac{-i \Omega^2}{a^2}\frac{1}{ \zeta} \frac{d V(\zeta)}{d \zeta} - \frac{i \Omega^2 V(\zeta)}{a^2} &= \frac{p_0}{\mu}\end{aligned}$$
#'If we divide by $\frac{-i \Omega_n^2}{a^2}$, then separate the solution
#'into homogenous and particular solutions: $$\begin{aligned}
#'    V &= V_h + V_p\\
#'    0 &= \frac{d^2 V(\zeta)}{d \zeta^2} + \frac{1}{\zeta}\frac{d V(\zeta)}{d \zeta} +  V(\zeta) \ \textrm{($V_h$)}\\
#'    0 &= \zeta^2 \frac{d^2 V(\zeta)}{d \zeta^2} + \zeta\frac{d V(\zeta)}{d \zeta} +  (\zeta^2 - \alpha^2) V(\zeta)  \textrm{, Where $\alpha = 0$}\\\end{aligned}$$
#'This last equation is known as Bessel's differential equation, and has a
#'solution of the form: $V_h(\zeta) = A_n J_0(\zeta) + B_n Y_0(\zeta)$,
#'where $J_0$ is the 0th order Bessel function of the 1st kind, and $Y_0$
#'is the 0th order Bessel function of the second kind. However,
#'$Y_0(0) = -\infty$, so for the solution to be finite at the center,
#'$B_n = 0$. $$\begin{aligned}
#'    V(\zeta) &= V_h(\zeta) + V_p(\zeta)\\
#'         &= A_n J_0(\zeta) + \frac{i p_o a^2}{\mu \Omega^2} \end{aligned}$$
#'Our next boundary condition is that at the wall, the velocity must be 0.
#'At the wall, $\zeta = \Lambda_n \frac{a}{a}$ $V(\Lambda_n) = 0$. If we
#'plug in this boundary condition and then substitute $V(\zeta)$ back into
#'$v_z = V(\zeta) e^{i n \omega t}$, we get: $$\begin{aligned}
#'    0 &= A_n J_0(\Lambda_n) + \frac{i p_o a^2}{\mu \Omega_n^2} \\
#'    A_n &= -\frac{\frac{i p_o a^2}{\mu \Omega_n^2}}{J_0(\Lambda_n)} \\
#'    V(\zeta) &= \frac{i p_o a^2}{\mu \Omega_n^2}\left(1 - \frac{J_0(\zeta)}{J_0(\Lambda_n)} \right)\\
#'    v_z^n(\zeta,t) &= \frac{i p_o a^2}{\mu \Omega_n^2}\left(1 - \frac{J_0(\zeta)}{J_0(\Lambda_n)} \right) e^{i n \omega t}\end{aligned}$$
#'
#'#### Flow Rate
#'
#'$$\begin{aligned}
#'    Q_n(t) &= \int_0^a 2 \pi r v_{z}^n(r,t) dr\\
#'           &= \frac{2 \pi p_o a^2}{\mu \Omega_n^2} e^{i n \omega t} \int_0^{a} r \left(1 - \frac{J_0(\zeta)}{J_o(\Lambda)} \right) dr\end{aligned}$$
#'We can substitute $r = \frac{a}{\Lambda} \zeta$ and
#'$dr = \frac{a}{\Lambda} d \zeta$. We can solve this integral by using
#'the identity: $\int_0^a x J_0(x) dx = a J_1(a)$ $$\begin{aligned}
#'           &= \frac{2 \pi p_o a^2}{\mu \Omega_n^2} e^{i n \omega t} \frac{a^2}{\Lambda_n^2 J_0(\Lambda_n)} \int_0^{\Lambda_n} \zeta \left(J_0(\Lambda_n) - J_0(\zeta) \right) d \zeta\\
#'           &= \frac{2 \pi p_o a^2}{\mu \Omega_n^2} e^{i n \omega t} \frac{a^2}{\Lambda_n^2 J_0(\Lambda_n)} \left(\int_0^{\Lambda_n} \zeta J_0(\Lambda_n) - \int_0^{\Lambda_n}\zeta J_0(\zeta) d \zeta \right)\\
#'           &= \frac{2 \pi p_o a^2}{\mu \Omega_n^2} e^{i n \omega t} \frac{a^2}{\Lambda_n^2 J_0(\Lambda_n)} \left(\frac{\Lambda_n^2}{2} J_0(\Lambda_n) - \Lambda_n J_1(\Lambda_n) \right)\\
#'           &= \frac{2 \pi p_o a^2}{\mu \Omega_n^2} e^{i n \omega t} \frac{a^2}{2} \left(1 - \frac{2 J_1(\Lambda_n)}{\Lambda_n J_0(\Lambda_n)}\right)\\
#'           &= \frac{\pi p_o a^4}{\mu \Omega_n^2} \left(1 - \frac{2 J_1(\Lambda_n)}{\Lambda_n J_0(\Lambda_n)}\right)e^{i n \omega t} \end{aligned}$$
#'
#'#### Shear Stress
#'
#'We can use the constitutive relation
#'$\tau_{rz} = 2 \mu \frac{1}{2}\left(    \frac{\partial v_z}{\partial r}+     \frac{\partial v_r}{\partial z}\right)$,
#'and the identity $\frac{d J_0(x)}{d x} = -J_1(x)$: $$\begin{aligned}
#'    \tau_n(t) &= -\mu \left(    \frac{\partial v_n(r,t)}{\partial r}\right)_{r=a} e^{i n \omega t}\\
#'              &= -\frac{i p_o a^2}{\Omega_n^2}\left[\frac{d}{d r}\left(1 - \frac{J_0(\zeta)}{J_0(\Lambda_n)} \right)   \right]_{r=a} e^{i n \omega t}\\
#'              &= -\frac{i p_o a^2}{\Omega_n^2}\frac{\Lambda_n}{a}\left[\frac{d}{d \zeta}\left(1 - \frac{J_0(\zeta)}{J_0(\Lambda_n)} \right)   \right]_{\zeta=\Lambda_n} e^{i n \omega t}\\
#'              &= -\frac{i p_o a^2}{\Omega_n^2}\frac{\Lambda_n}{a}\left[\left(\frac{J_1(\zeta)}{J_0(\Lambda_n)} \right)   \right]_{\zeta=\Lambda_n} e^{i n \omega t}\\
#'              &= -\frac{i p_o a^2}{\Omega_n^2}\frac{\Lambda_n}{a}\left(\frac{J_1(\Lambda_n)}{J_0(\Lambda_n)} \right) e^{i n \omega t}\\
#'              &= -\frac{p_o a}{\Lambda_n}\left(\frac{J_1(\Lambda_n)}{J_0(\Lambda_n)} \right) e^{i n \omega t}\\\end{aligned}$$
#'The last simplification was made since
#'$\frac{i \Lambda}{\Omega^2} = \frac{1}{\Lambda}$
#'
#'Pulsatile
#'=========
#'
#'#### Velocity
#'
#'Flow is driven by the following combined pressure gradient:
#'$$\frac{d p}{d z} = -p_s - p_o \sum_{n=1}^N C_n e^{i (n \omega t + \phi_n)} \textrm{, Where $\omega = \frac{2 \pi}{T}$}$$
#'Combining our steady and oscillatory velocity profiles, we get:
#'$$v_{z,o}(r,t) = \frac{p_s a^2}{4 \mu} \left(1 - \left(\frac{r}{a}\right)^2\right) + p_o \displaystyle \sum_{n = 1}^{N} \frac{i C_n a^2}{\mu \Omega_n^2} \left(1 - \frac{J_0 (\zeta)}{J_0 (\Lambda_n)}\right) e^{i (\omega n t - \phi_n)}$$
#'
#'#### Flow Rate
#'
#'Similarly, combining the steady and oscillatory portions of flow rate,
#'we get:
#'$$Q(t) = \frac{\pi a^4 p_s}{8 \mu} +\displaystyle \sum_{n=1}^N\frac{\pi p_n a^4}{\mu \Omega_n^2} \left(1 - \frac{2 J_1(\Lambda_n)}{\Lambda_n J_0(\Lambda_n)}\right)e^{i n \omega t}$$
#'
#'#### Wall shear stress
#'
#'And finally, for wall shear stress we can similarly combine our steady
#'and oscillatory profiles.
#'$$\tau(t) = -\frac{a p_s}{2} +  \sum_{n=1}^{N}-\frac{p_n a}{\Lambda_n}\left(\frac{J_1(\Lambda_n)}{J_0(\Lambda_n)} \right) e^{i n \omega t}\\
#'$$
# Don't show default figure, only html5 video
#+ results='raw'
#!/usr/bin/env python
import numpy as np
import re
from scipy.special import jv
import matplotlib.animation as animation
import matplotlib.pyplot as plt
from IPython.display import display_html
def womersley(rho,omega,mu,a): # Oscillatory portion
    return a*np.sqrt((rho*omega)/mu)
def reynolds(rho, u_max, D, mu): # steady portion
    return (rho*u_max_D)/mu

def pressure(dp_dz,omega,waves,t):
    """TODO: Docstring for pressure.

    Parameters
    ----------
    dp_dz   : TODO
    omega   : TODO
    waves   : TODO
    t       : TODO

    Returns
    -------
    TODO

    """
    c_ns,phi_ns,omegas,lambdas = waves
    press_ind = np.array([c_n*np.exp(1j*(omega*(n+1)*t + phi_n))
            for n,(c_n,phi_n) in enumerate(zip(c_ns,phi_ns))])
    press = np.sum(press_ind,axis=0) + dp_dz
    return press.real
def Q(a,omega,mu,waves,dp_dz,t):
    """TODO: Docstring for Q.

    Parameters
    ----------
    a       : TODO
    omega   : TODO
    mu      : TODO
    waves   : TODO
    dp_dz   : TODO
    t       : TODO

    Returns
    -------
    TODO

    """
    q_steady = (dp_dz*np.pi*a**4)/(8*mu)
    c_ns,phi_ns,omegas,lambdas = waves
    omegas = np.array([np.sqrt((n*rho*omega)/mu)*a for n in range(1,len(c_ns)+1)])
    lambdas = ((1j - 1)/np.sqrt(2))*omegas
    q_unsteady = np.array([((1j*np.pi*c_n*a**4)/(mu*om**2))*\
            (1 - (2*jv(1,lam))/(lam*jv(1,lam)))*\
            np.exp(1j*((n + 1)*omega*t + phi_n))\
            for n,(c_n,phi_n,om,lam) in enumerate(zip(c_ns,phi_ns,omegas,lambdas))])
    q_total = np.sum(q_unsteady,axis=0) + q_steady
    return q_total.real
def tau(a,waves,dp_dz,t):
    """TODO: Docstring for tau.

    Parameters
    ----------
    a       : TODO
    waves   : TODO
    dp_dz   : TODO
    t       : TODO

    Returns
    -------
    TODO

    """
    tau_steady = dp_dz*a/2
    c_ns,phi_ns,omegas,lambdas = waves
    tau_unsteady = np.array([(-c_n*a/lam)*(jv(1,lam)/jv(0,lam))*\
            np.exp(1j*((n + 1)*omega*t + phi_n))\
            for n,[c_n,phi_n,lam] in enumerate(zip(c_ns,phi_ns,lambdas))])
    tau_tot = np.sum(tau_unsteady,axis=0) + tau_steady
    return tau_tot.real
def flow(a,omega,mu,rho,waves,dp_dz,animated):
    """TODO: Docstring for flow.
        Parameters
        ----------
        a           : float
            Radius of the inside of the tube
        omega       : float
            Oscillation frequency in radians/s
        mu          : float
            Newtonian viscosity of fluid
        rho         : float
            Constant density of fluid
        waves       : tuple of the form (list(C_n),list(Phi_n))
            2-Tuple of lists of the parameters C_n and Phi_n
        dp_dz       : float
            Steady state pressure grad
        animated    : boolean
            Show animated movie or phased graphs.

        Returns
        -------
        TODO
    """
    plt.style.use('ggplot')
    fig = plt.figure(figsize=(10,15))
    ax_v = fig.add_subplot(411)
    ax_p = fig.add_subplot(412)
    ax_q = fig.add_subplot(413,sharex=ax_p)
    ax_t = fig.add_subplot(414,sharex=ax_p)

    #Plot "tube" on velocity profile graph
    ax_v.plot([-250,250],[1,1])
    ax_v.plot([-250,250],[-1,-1])
    ax_v.set_xlim([-200,200])

    #TODO: set x_limit for all other axes
    # 2 cycles
    times = np.linspace(0,2*T,1000)

    # Set Omega, Lambda, Zeta, C_n values, Phi_n values, and the constant coefficient out front
    r_a = np.linspace(-1,1)
    c_ns,phi_ns,omegas,lambdas = waves

    # Steady state flow rate
    v_steady = dp_dz*a**2/(4*mu)*(1 - r_a**2)

    # Pressure
    flow_pressure = pressure(dp_dz,omega,waves,times)
    flow_Q = Q(a,omega,mu,waves,dp_dz,times)
    flow_tau = tau(a,waves,dp_dz,times)
    ax_p.plot(times,flow_pressure)
    ax_q.plot(times,flow_Q)
    ax_t.plot(times,flow_tau)

    # X and Y Labels
    ax_v.set_xlabel(r'$v_z(r/a)$')
    ax_t.set_xlabel('time (s)')

    ax_v.set_ylabel(r'Normalized Radius ($r/a$)')
    ax_p.set_ylabel(r'Pressure ($dyne/cm^2$)')
    ax_q.set_ylabel(r'Flow Rate ($cm^3/s$)')
    ax_t.set_ylabel(r'Wall shear stress ($dyne/cm^2$)')

    #Limits
    ax_p.set_xlim([0,times[-1]])
    p_range = (np.max(flow_pressure)-np.min(flow_pressure))/2
    q_range = (np.max(flow_Q)-np.min(flow_Q))/2
    t_range = (np.max(flow_tau)-np.min(flow_tau))/2
    ax_p.set_ylim([np.min(flow_pressure)-p_range,np.max(flow_pressure)+p_range])
    ax_q.set_ylim([np.min(flow_Q)-q_range,np.max(flow_Q)+q_range])
    ax_t.set_ylim([np.min(flow_tau)-t_range,np.max(flow_tau)+t_range])

    #Base lines:
    p_follow, = ax_p.plot([0,0],[np.min(flow_pressure)-p_range*1.1,np.max(flow_pressure)+p_range*1.1])
    ax_p.plot([-10,20],[50,50])
    q_follow, = ax_q.plot([0,0],[np.min(flow_Q)-q_range*1.1,np.max(flow_Q)+q_range*1.1])
    t_follow, = ax_t.plot([0,0],[np.min(flow_tau)-t_range*1.1,np.max(flow_tau)+t_range*1.1])
    u_an = np.array([((1j*c_n*a**2)/(mu*om**2))*(1 - (jv(0,abs(r_a)*lam)/jv(0,lam)))* \
            np.exp(1j*((n + 1)*omega*0 + phi_n)) \
            for n,[c_n,phi_n,om,lam] in enumerate(zip(c_ns,phi_ns,omegas,lambdas))])
    u_a_tot = np.sum(u_an,axis=0) #TODO + u_steady
    line, = ax_v.plot(u_a_tot.real,r_a)

    def animate(t):
        v_unsteady = np.array([((1j*c_n*a**2)/(mu*om**2))*(1 - (jv(0,abs(r_a)*lam)/jv(0,lam)))* \
                np.exp(1j*((n + 1)*omega*t + phi_n)) \
                for n,[c_n,phi_n,om,lam] in enumerate(zip(c_ns,phi_ns,omegas,lambdas))])
        v_pulse = np.sum(v_unsteady,axis=0) + v_steady
        p_follow.set_xdata([t, t])
        q_follow.set_xdata([t, t])
        t_follow.set_xdata([t, t])

        line.set_xdata(v_pulse.real)  # update the data.
        return line,
    def frames_gen():
        # 100 frames from 0 to 2 periods.
        T = 1/(omega/(2*np.pi))
        for frame in np.linspace(0,2*T,100):
            yield frame
    # 100 frames from 0 to 2 periods.  100 frames with 75ms per frame we have a 7.5s long animation with a 1s delay at the end
    ani = animation.FuncAnimation(
            fig, animate, interval=75, blit=False, frames=frames_gen,repeat_delay=1000)
    vid = ani.to_html5_video()
    # Change width to 100% and remove height
    vid = re.sub(r'width="\d+"',r'width="100%"',vid)
    vid = re.sub(r'height="\d+"',r'',vid)
    #display video
    display_html(vid,raw=True)
#' We can first set our constants, namely $\rho = 1.06\ g/cm^3$, $\mu = 4 \cP$, $C_n$ values,  $\Phi_n$ values, and $T = 0.8$
T = .8
omega = 2*np.pi/T #1.25 hertz
rho = 1.06 #1060 kg/m^2 == 1.06g/cm^3
mu = .04
c_n = np.array([7.58,5.41,1.52,.52,.83,.69,.26,.54,.27,.1])
phi_n = np.array([-174,89,-22,-34,-127,135,152,44,-72,11])*np.pi/180

#' Now we visualize pulsatile flow with the first set of results with $a = .65 \ cm$, $C = 50 \ dynes/cm^3$, and $D = 10\ dynes/cm^3$
#+ results='raw',fig = False
a = .65 #.65 cm
omegas = np.array([np.sqrt((n*rho*omega)/mu)*a for n in range(1,len(c_n)+1)])
lambdas = ((1j - 1)/np.sqrt(2))*omegas
waves = (10*c_n,phi_n,omegas,lambdas)
dp_dz = 50
flow(a,omega,mu,rho,waves,dp_dz,True)

#' Now we visualize pulsatile flow with the first set of results with $a = .3 \ cm$, $C = 20 \ dynes/cm^3$, and $D = 30\ dynes/cm^3$
#+ results='raw',fig = False
a = .3 #.65 cm
omegas = np.array([np.sqrt((n*rho*omega)/mu)*a for n in range(1,len(c_n)+1)])
lambdas = ((1j - 1)/np.sqrt(2))*omegas
waves = (30*c_n,phi_n,omegas,lambdas)
dp_dz = 20
flow(a,omega,mu,rho,waves,dp_dz,True)
