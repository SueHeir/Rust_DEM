use nalgebra::Vector3;

pub struct Box {
    pub position: Vector3<i32>,
    pub real: Vec<i32>,
    pub ghost: Vec<i32>,
    pub lo: Vector3<f64>,
    pub hi: Vector3<f64>,
}

impl Box {
    pub fn is_position_in_box(&self, position: Vector3<f64>) -> bool {
        if position[0] >= self.lo[0]
            && position[0] < self.hi[0]
            && position[1] >= self.lo[1]
            && position[1] < self.hi[1]
            && position[2] >= self.lo[2]
            && position[2] < self.hi[2]
        {
            return true;
        }

        return false;
    }

    pub fn is_position_in_max_radius_enlarged_box(
        &self,
        position: Vector3<f64>,
        max_radius: f64,
    ) -> bool {
        if position[0] >= self.lo[0] - max_radius
            && position[0] <= self.hi[0] + max_radius
            && position[1] >= self.lo[1] - max_radius
            && position[1] <= self.hi[1] + max_radius
            && position[2] >= self.lo[2] - max_radius
            && position[2] <= self.hi[2] + max_radius
        {
            return true;
        }

        return false;
    }

    pub fn _is_sphere_aabb_in_box(&self, position: Vector3<f64>, radius: f64) -> bool {
        if self.is_position_in_box(position + Vector3::new(radius, radius, radius)) {
            return true;
        } else if self.is_position_in_box(position + Vector3::new(radius, radius, -radius)) {
            return true;
        } else if self.is_position_in_box(position + Vector3::new(radius, -radius, radius)) {
            return true;
        } else if self.is_position_in_box(position + Vector3::new(radius, -radius, -radius)) {
            return true;
        } else if self.is_position_in_box(position + Vector3::new(-radius, radius, radius)) {
            return true;
        } else if self.is_position_in_box(position + Vector3::new(-radius, radius, -radius)) {
            return true;
        } else if self.is_position_in_box(position + Vector3::new(-radius, -radius, radius)) {
            return true;
        } else if self.is_position_in_box(position + Vector3::new(-radius, -radius, -radius)) {
            return true;
        }

        return false;
    }

    pub fn is_sphere_aabb_in_radius_enlarged_box(
        &self,
        position: Vector3<f64>,
        radius: f64,
        max_radius: f64,
    ) -> bool {
        if self.is_position_in_max_radius_enlarged_box(
            position + Vector3::new(radius, radius, radius),
            max_radius,
        ) {
            return true;
        } else if self.is_position_in_max_radius_enlarged_box(
            position + Vector3::new(radius, radius, -radius),
            max_radius,
        ) {
            return true;
        } else if self.is_position_in_max_radius_enlarged_box(
            position + Vector3::new(radius, -radius, radius),
            max_radius,
        ) {
            return true;
        } else if self.is_position_in_max_radius_enlarged_box(
            position + Vector3::new(radius, -radius, -radius),
            max_radius,
        ) {
            return true;
        } else if self.is_position_in_max_radius_enlarged_box(
            position + Vector3::new(-radius, radius, radius),
            max_radius,
        ) {
            return true;
        } else if self.is_position_in_max_radius_enlarged_box(
            position + Vector3::new(-radius, radius, -radius),
            max_radius,
        ) {
            return true;
        } else if self.is_position_in_max_radius_enlarged_box(
            position + Vector3::new(-radius, -radius, radius),
            max_radius,
        ) {
            return true;
        } else if self.is_position_in_max_radius_enlarged_box(
            position + Vector3::new(-radius, -radius, -radius),
            max_radius,
        ) {
            return true;
        }

        return false;
    }

    pub fn is_periodic_sphere(
        &self,
        position: Vector3<f64>,
        radius: f64,
        max_radius: f64,
        d_data: &DomainData,
    ) -> bool {
        if position[0] - radius <= 0.0 + max_radius {
            return true;
        }
        if position[0] + radius >= d_data.domain[0] - max_radius {
            return true;
        }
        if position[1] - radius <= 0.0 + max_radius {
            return true;
        }
        if position[1] + radius >= d_data.domain[1] - max_radius {
            return true;
        }
        if position[2] - radius <= 0.0 + max_radius {
            return true;
        }
        if position[2] + radius >= d_data.domain[2] - max_radius {
            return true;
        }

        return false;
    }
}

#[derive()]
pub struct DomainData {
    pub(crate) domain: Vector3<f64>,
    pub(crate) domain_volume: f64,
    pub(crate) collision_boxes: Vector3<i32>,
    pub(crate) lees_edwards_boundary: f64,
    pub(crate) g_data: Vec<Vec<Vec<Box>>>,
    pub(crate) gravity: Vector3<f64>,
}
