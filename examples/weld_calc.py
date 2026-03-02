import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

class Weld():
    """
    Class representing a single weld bead, now supporting any orientation.
    """
    _next_id = 1  # class counter for default names

    def __init__(self, name=None):
        self.h = 0.005                     # weld throat size (m)
        self.start_point = np.array([0, 0, 0], dtype=float)
        self.end_point = np.array([1, 0, 0], dtype=float)
        self.l = None                      # length (will be computed)
        self.color = "blue"
        self.line_width = 3
        # Assign name
        if name is None:
            self.name = f"Weld {Weld._next_id}"
        else:
            self.name = name
        Weld._next_id += 1

    def weld_calcs(self):
        """
        Compute geometric properties of the weld, including local moments of inertia
        rotated to global axes.
        """
        # Length and direction vector
        self.l = np.linalg.norm(self.start_point - self.end_point)
        direction = (self.end_point - self.start_point) / self.l
        # Angle with respect to global X axis (in radians)
        self.theta = np.arctan2(direction[1], direction[0])

        # Area of throat (rectangle of dimensions l x (h/√2))
        throat_width = self.h / np.sqrt(2.)
        self.A = throat_width * self.l
        self.centroid = (self.start_point + self.end_point) * 0.5

        # Local moments of inertia (about centroid of the throat rectangle)
        # Axis along the weld: I_long = (1/12) * throat_width * l^3
        # Axis perpendicular to weld: I_trans = (1/12) * l * throat_width^3
        I_long = (1./12.) * throat_width * self.l**3
        I_trans = (1./12.) * self.l * throat_width**3

        # Rotate to global axes using transformation formulas
        c = np.cos(self.theta)
        s = np.sin(self.theta)
        # I_ξ = I_long, I_η = I_trans
        self.Ixx_loc = I_long * s**2 + I_trans * c**2
        self.Iyy_loc = I_long * c**2 + I_trans * s**2
        self.Ixy_loc = (I_long - I_trans) * s * c

        # JG (polar moment) remains as before (sum of I_long + I_trans)
        self.JG = I_long + I_trans

        # For backward compatibility
        self.IxG = self.Ixx_loc
        self.IyG = self.Iyy_loc

    def print_report(self):
        """Print a summary of the weld properties."""
        if self.l is None:
            print("Run the weld_calcs method first")
        else:
            print(f"""############
--- REPORT ---
############
Name  = {self.name}
l     = {self.l:.6f} m
Area  = {self.A:.6f} m²
JG    = {self.JG:.6f} m⁴
Ixx_loc = {self.Ixx_loc:.6f} m⁴
Iyy_loc = {self.Iyy_loc:.6f} m⁴
Ixy_loc = {self.Ixy_loc:.6f} m⁴
Centroid:
xG    = {self.centroid[0]:.6f} m
yG    = {self.centroid[1]:.6f} m
Angle = {np.degrees(self.theta):.2f}°""")

    def plot_weld(self):
        """Plot the individual weld bead."""
        plt.figure(figsize=(3, 3))
        plt.plot(np.vstack((self.start_point, self.end_point))[:,0],
                 np.vstack((self.start_point, self.end_point))[:,1],
                 linewidth=self.line_width, color=self.color)
        plt.grid(True)
        plt.axis('equal')
        plt.show()

    def __str__(self):
        return "It is a class called Weld, represented by two points and a width"


class Weld_joint():
    """
    Class representing a joint composed of several welds.
    Handles welds of any orientation using full inertia tensor.
    """
    def __init__(self):
        self.welds = None
        self.figure_size = (4, 4)
        self.figure_name = "my_weld_joint.png"
        self.dpi = 100

    def translate_forces_to_centroid(self, forces, points_F):
        """
        Translate multiple forces to the joint centroid.
        """
        F_total = np.zeros(3)
        M_total = np.zeros(3)

        for i, F in enumerate(forces):
            F_total += F
            r = points_F[i] - self.centroid
            M_total += np.cross(r, F)

        return F_total, M_total

    def weld_calcs(self):
        """
        Compute geometric properties of the joint:
        - Total area
        - Centroid
        - JG (polar moment for torsion)
        - Ixx, Iyy, Ixy (full inertia tensor for bending)
        """
        n = len(self.welds)
        self.r = np.ones(n)
        self.A = 0.
        self.centroid = np.array([0, 0, 0], dtype=float)
        self.JG = 0.
        self.Ixx = 0.
        self.Iyy = 0.
        self.Ixy = 0.

        # First pass: compute centroid
        for weld in self.welds:
            self.A += weld.A
            self.centroid[0] += weld.A * weld.centroid[0]
            self.centroid[1] += weld.A * weld.centroid[1]

        self.centroid[0] *= 1. / self.A
        self.centroid[1] *= 1. / self.A

        # Second pass: accumulate JG, Ixx, Iyy, Ixy using Steiner
        for i, weld in enumerate(self.welds):
            dx = weld.centroid[0] - self.centroid[0]
            dy = weld.centroid[1] - self.centroid[1]
            r = np.sqrt(dx**2 + dy**2)
            self.r[i] = r

            self.JG += weld.JG + weld.A * r**2
            self.Ixx += weld.Ixx_loc + weld.A * dy**2
            self.Iyy += weld.Iyy_loc + weld.A * dx**2
            self.Ixy += weld.Ixy_loc + weld.A * dx * dy

    def weld_evaluation(self, forces, points_F, weld_points):
        """
        Evaluate stresses at given points on the welds for multiple loads.
        Uses full bending formula for normal stresses.
        """
        F_total, M_total = self.translate_forces_to_centroid(forces, points_F)
        # Simply call the centroid-based version
        self.weld_evaluation_from_centroid_loads(F_total, M_total, weld_points)

    # ========= NUEVO MÉTODO =========
    def weld_evaluation_from_centroid_loads(self, F_total, M_total, weld_points):
        """
        Evaluate stresses at given points on the welds for loads already
        reduced to the centroid.

        Parameters:
        -----------
        F_total : array (3,)
            Resultant force vector at the centroid [Fx, Fy, Fz] (N).
        M_total : array (3,)
            Resultant moment vector about the centroid [Mx, My, Mz] (N·m).
        weld_points : array (n x 3)
            Points where stresses are evaluated.
        """
        self.F_total = F_total
        self.M_total = M_total

        Fx, Fy, Fz = F_total
        Mx, My, Mz = M_total

        n_points = len(weld_points)
        self.TAU_1 = np.zeros((n_points, 3))   # primary shear
        self.TAU_2 = np.zeros((n_points, 3))   # torsional shear (only Mz)
        self.TAU = np.zeros((n_points, 3))     # total shear
        self.SIGMA = np.zeros(n_points)         # normal stress
        self.VON_MISES = np.zeros(n_points)     # von Mises equivalent stress

        M_torsion = np.array([0, 0, Mz])        # use only Mz for torsion

        # Precompute denominator for bending formula (det = Ixx*Iyy - Ixy^2)
        det = self.Ixx * self.Iyy - self.Ixy**2
        if abs(det) < 1e-12:
            use_simple_bending = True
        else:
            use_simple_bending = False

        for i, point in enumerate(weld_points):
            # Primary shear (Fx, Fy)
            self.TAU_1[i] = -np.array([Fx, Fy, 0]) / self.A

            # Torsional shear (Mz)
            r = point - self.centroid
            self.TAU_2[i] = np.cross(-M_torsion, r) / self.JG

            self.TAU[i] = self.TAU_1[i] + self.TAU_2[i]

            # Normal stress (due to Fz, Mx, My)
            sigma_axial = Fz / self.A if self.A != 0 else 0

            if use_simple_bending:
                sigma_bx = Mx * (point[1] - self.centroid[1]) / self.Ixx if self.Ixx != 0 else 0
                sigma_by = My * (point[0] - self.centroid[0]) / self.Iyy if self.Iyy != 0 else 0
                sigma_bending = sigma_bx + sigma_by
            else:
                x = point[0] - self.centroid[0]
                y = point[1] - self.centroid[1]
                sigma_bending = ( (Mx * self.Iyy - My * self.Ixy) * y
                                 + (My * self.Ixx - Mx * self.Ixy) * x ) / det

            self.SIGMA[i] = sigma_axial + sigma_bending

            # von Mises equivalent stress
            tau_mag = np.linalg.norm(self.TAU[i])
            self.VON_MISES[i] = np.sqrt(self.SIGMA[i]**2 + 3 * tau_mag**2)

    def weld_point_evaluation(self, forces, points_F, point):
        """
        Evaluate stresses at a single point (e.g., a critical corner).
        """
        point_array = np.array([point])
        self.weld_evaluation(forces, points_F, point_array)

        self.tau_1 = self.TAU_1[0]
        self.tau_2 = self.TAU_2[0]
        self.tau = self.TAU[0]
        self.sigma = self.SIGMA[0]
        self.von_mises = self.VON_MISES[0]

    def plot_welds(self):
        """Plot all welds and the joint centroid."""
        plt.figure(figsize=self.figure_size)
        for weld in self.welds:
            plt.plot(np.vstack((weld.start_point, weld.end_point))[:,0],
                     np.vstack((weld.start_point, weld.end_point))[:,1],
                     linewidth=weld.line_width, color=weld.color)

        plt.plot(self.centroid[0], self.centroid[1],
                 marker="o", markeredgecolor="red",
                 markersize=6, fillstyle='none', markeredgewidth=3)
        plt.grid(True)
        plt.axis('equal')
        plt.show()

    def plot_vectors(self, weld_points, plot_type='tau', scale_vectors=1,
                     scale_magnitud=1e-6, scale_points=1000, show_grid=True,
                     save_figure=False, plot_TAU_1=False, plot_TAU_2=False,
                     show_weld_lines=True, cmap="jet"):
        """
        Visualise stress results with vectors and colour maps.
        """
        X = weld_points[:, 0] * scale_points
        Y = weld_points[:, 1] * scale_points

        plt.figure(figsize=self.figure_size)

        if show_weld_lines:
            for weld in self.welds:
                plt.plot(np.vstack((weld.start_point, weld.end_point))[:,0] * scale_points,
                         np.vstack((weld.start_point, weld.end_point))[:,1] * scale_points,
                         linewidth=weld.line_width, color=weld.color, alpha=0.5,
                         zorder=1)

        plt.plot(self.centroid[0] * scale_points, self.centroid[1] * scale_points,
                 marker="o", markeredgecolor="red", markersize=6,
                 fillstyle='none', markeredgewidth=3, zorder=2)

        if plot_type == 'tau':
            magnitud = np.linalg.norm(self.TAU, axis=1) * scale_magnitud
            normalized = (magnitud - np.min(magnitud)) / (np.max(magnitud) - np.min(magnitud) + 1e-12)

            for i in range(len(self.TAU)):
                plt.arrow(X[i], Y[i],
                          self.TAU[i, 0] * scale_magnitud,
                          self.TAU[i, 1] * scale_magnitud,
                          head_width=3*scale_vectors, head_length=3*scale_vectors,
                          color=plt.cm.jet(normalized[i]),
                          zorder=3)

            norm = Normalize(vmin=np.min(magnitud), vmax=np.max(magnitud))
            mappable = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
            mappable.set_array([])
            plt.colorbar(mappable, ax=plt.gca(), label='τ (MPa)')

        elif plot_type == 'sigma':
            sigma_mpa = self.SIGMA * 1e-6
            scatter = plt.scatter(X, Y, c=sigma_mpa, cmap=cmap,
                                   s=50, edgecolors='black', linewidth=0.5,
                                   zorder=4)
            plt.colorbar(scatter, label='σ (MPa)')

        elif plot_type == 'von_mises':
            vm_mpa = self.VON_MISES * 1e-6
            scatter = plt.scatter(X, Y, c=vm_mpa, cmap=cmap,
                                   s=50, edgecolors='black', linewidth=0.5,
                                   zorder=4)
            plt.colorbar(scatter, label='σ_von Mises (MPa)')

        if plot_TAU_1 and hasattr(self, 'TAU_1'):
            for i in range(len(self.TAU_1)):
                plt.arrow(X[i], Y[i],
                          self.TAU_1[i, 0] * scale_magnitud,
                          self.TAU_1[i, 1] * scale_magnitud,
                          head_width=2*scale_vectors, head_length=2*scale_vectors,
                          linestyle='--', alpha=0.5, color='gray',
                          zorder=3)

        if plot_TAU_2 and hasattr(self, 'TAU_2'):
            for i in range(len(self.TAU_2)):
                plt.arrow(X[i], Y[i],
                          self.TAU_2[i, 0] * scale_magnitud,
                          self.TAU_2[i, 1] * scale_magnitud,
                          head_width=2*scale_vectors, head_length=2*scale_vectors,
                          linestyle=':', alpha=0.5, color='gray',
                          zorder=3)

        plt.xlabel('X (mm)')
        plt.ylabel('Y (mm)')
        plt.grid(show_grid)
        plt.axis('equal')

        if save_figure:
            plt.savefig(self.figure_name, dpi=self.dpi, bbox_inches='tight')
        else:
            plt.show()

    def safety_factor(self, yield_strength, criterion='vm'):
        """
        Calculate the safety factor at each evaluated point.
        """
        if not hasattr(self, 'VON_MISES'):
            raise ValueError("You must run weld_evaluation() first")

        if criterion.lower() in ('vm', 'vonmises', 'von mises'):
            equiv_stress = self.VON_MISES
        elif criterion.lower() in ('tresca', 'max_shear'):
            tau_mag = np.linalg.norm(self.TAU, axis=1)
            equiv_stress = 2 * np.sqrt((self.SIGMA/2)**2 + tau_mag**2)
        else:
            raise ValueError("Invalid criterion. Use 'vm' or 'tresca'")

        n = yield_strength / equiv_stress
        n_min = np.min(n)
        return n_min, n

    def plot_sigma_unfolded(self, weld_points, labels_cordones=None, figsize=(10,4)):
        """
        Bar plot of normal stress (σ) along the unfolded joint.
        """
        if not hasattr(self, 'SIGMA'):
            raise ValueError("You must run weld_evaluation() first")

        n_points = len(weld_points)
        weld_indices = np.zeros(n_points, dtype=int)

        # Assign each point to the nearest weld (2D distance to segment)
        for j, point in enumerate(weld_points):
            p = point[:2]
            min_dist = float('inf')
            best_i = -1
            for i, weld in enumerate(self.welds):
                A = weld.start_point[:2]
                B = weld.end_point[:2]
                AB = B - A
                L = weld.l
                if L == 0:
                    dist = np.linalg.norm(p - A)
                else:
                    t = np.dot(p - A, AB) / (L*L)
                    if t < 0:
                        closest = A
                    elif t > 1:
                        closest = B
                    else:
                        closest = A + t * AB
                    dist = np.linalg.norm(p - closest)
                if dist < min_dist:
                    min_dist = dist
                    best_i = i
            weld_indices[j] = best_i

        sorted_sigma = []
        bar_colors = []
        boundaries = [0]

        for i, weld in enumerate(self.welds):
            mask = weld_indices == i
            pts_i = weld_points[mask]
            if len(pts_i) == 0:
                boundaries.append(boundaries[-1])
                continue

            # Curvilinear coordinate s along the weld
            A = weld.start_point[:2]
            B = weld.end_point[:2]
            AB = B - A
            L = weld.l
            s_i = []
            for p in pts_i:
                p2d = p[:2]
                t = np.dot(p2d - A, AB) / (L*L) if L != 0 else 0
                s_i.append(t * L)

            # Sort by s
            idx_sort = np.argsort(s_i)
            sigma_i = self.SIGMA[mask][idx_sort] * 1e-6  # convert to MPa
            sorted_sigma.extend(sigma_i)
            bar_colors.extend([weld.color] * len(sigma_i))
            boundaries.append(boundaries[-1] + len(sigma_i))

        plt.figure(figsize=figsize)
        x_pos = np.arange(len(sorted_sigma))
        plt.bar(x_pos, sorted_sigma, color=bar_colors, alpha=0.7)

        if labels_cordones is None:
            labels_cordones = [weld.name for weld in self.welds]

        for idx, (start, end) in enumerate(zip(boundaries[:-1], boundaries[1:])):
            if end - start > 0:
                plt.axvline(end - 0.5, color='gray', linestyle='--', alpha=0.5)
                center = (start + end - 1) / 2
                plt.text(center, plt.ylim()[1]*0.9, labels_cordones[idx],
                         ha='center', va='center', bbox=dict(facecolor='white', alpha=0.8))

        plt.ylabel('Normal stress σ (MPa)')
        plt.xlabel('Points along the joint')
        plt.grid(axis='y', linestyle=':', alpha=0.6)
        plt.tight_layout()
        plt.show()

    def plot_n_unfolded(self, weld_points, yield_strength, criterion = 'vm',
                        labels_cordones=None, figsize=(10,4)):
        """
        Bar plot of safety factor (n) along the unfolded joint.
        """
        if not hasattr(self, 'VON_MISES'):
            raise ValueError("You must run weld_evaluation() first")

        _, n_values = self.safety_factor(yield_strength, criterion)

        n_points = len(weld_points)
        weld_indices = np.zeros(n_points, dtype=int)

        for j, point in enumerate(weld_points):
            p = point[:2]
            min_dist = float('inf')
            best_i = -1
            for i, weld in enumerate(self.welds):
                A = weld.start_point[:2]
                B = weld.end_point[:2]
                AB = B - A
                L = weld.l
                if L == 0:
                    dist = np.linalg.norm(p - A)
                else:
                    t = np.dot(p - A, AB) / (L*L)
                    if t < 0:
                        closest = A
                    elif t > 1:
                        closest = B
                    else:
                        closest = A + t * AB
                    dist = np.linalg.norm(p - closest)
                if dist < min_dist:
                    min_dist = dist
                    best_i = i
            weld_indices[j] = best_i

        sorted_n = []
        bar_colors = []
        boundaries = [0]

        for i, weld in enumerate(self.welds):
            mask = weld_indices == i
            pts_i = weld_points[mask]
            if len(pts_i) == 0:
                boundaries.append(boundaries[-1])
                continue

            A = weld.start_point[:2]
            B = weld.end_point[:2]
            AB = B - A
            L = weld.l
            s_i = []
            for p in pts_i:
                p2d = p[:2]
                t = np.dot(p2d - A, AB) / (L*L) if L != 0 else 0
                s_i.append(t * L)

            idx_sort = np.argsort(s_i)
            n_i = n_values[mask][idx_sort]
            sorted_n.extend(n_i)
            bar_colors.extend([weld.color] * len(n_i))
            boundaries.append(boundaries[-1] + len(n_i))

        plt.figure(figsize=figsize)
        x_pos = np.arange(len(sorted_n))
        plt.bar(x_pos, sorted_n, color=bar_colors, alpha=0.7)

        if labels_cordones is None:
            labels_cordones = [weld.name for weld in self.welds]

        for idx, (start, end) in enumerate(zip(boundaries[:-1], boundaries[1:])):
            if end - start > 0:
                plt.axvline(end - 0.5, color='gray', linestyle='--', alpha=0.5)
                center = (start + end - 1) / 2
                plt.text(center, plt.ylim()[1]*0.9, labels_cordones[idx],
                         ha='center', va='center', bbox=dict(facecolor='white', alpha=0.8))

        plt.axhline(y=1, color='red', linestyle='-', linewidth=1.5, label='n = 1')
        plt.xlabel('Points along the joint')
        plt.ylabel('SF')
        plt.grid(axis='y', linestyle=':', alpha=0.6)
        plt.legend()
        plt.tight_layout()
        plt.show()


# Auxiliary functions
def equidistant_points(weld, n):
    """Generate n equidistant points along a weld."""
    start = weld.start_point
    end = weld.end_point
    total_dist = np.linalg.norm(end - start)
    if total_dist == 0:
        return np.tile(start, (n, 1))
    direction = (end - start) / total_dist
    step = total_dist / (n - 1)
    points = np.zeros((n, 3))
    for i in range(n):
        points[i] = start + i * step * direction
    return points

def remove_repeats(points_matrix):
    """Remove duplicate points from a matrix."""
    unique = set(map(tuple, points_matrix))
    return np.array(list(unique))