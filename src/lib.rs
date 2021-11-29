use morpho::{*};

// This looks strange, but Particles can be quite arbitrary, for instance they can contain
// only `x' in one implementation or only `vx' in another.
trait Particle {}

// Something like radiative losses caused by the photon emission
trait SimplePusher<P: Particle> {
    fn push(&self, particle: P) -> P;
}

// Something like the Vay pusher, in that photon emission can be embedded in the middle point
trait CompositePusher<P: Particle>: SimplePusher<P> {
    fn with(&self, pusher: SimplePusher<P>) -> SimplePusher<P>;
}

/*trait Particle3r: Particle {
    fn r(&self) -> &V3;
}*/

// trait for particles...
trait RelativisticParticle3p: Particle {
    fn p(&self) -> &Vector3;
    fn g(&self) -> f64; // or just &f64 is better here?
                        // for cache f64 should be better?
                        // и по смыслу g просто вычисляется через p
}


/* trait for particles...
trait RelativisticParticle3r3p<V3: Vector3>: Particle {
    fn r(&self) -> &V3;
    fn p(&self) -> &V3;
    fn g(&self) -> &f64;
}*/

struct SynchrotronPusher {
    pub chi: f64,
}

impl<P: RelativisticParticle3p> SimplePusher<P> for SynchrotronPusher {
    fn push(&self, particle: P) -> P {
        return particle;
    }
}

// Implementation. Now we define what the particles are in the memory

/*struct RelativisticParticle_3r3p { r: Vector_3, p: Vector_3, g: f64 }

impl Particle for RelativisticParticle_3r3p {}

impl RelativisticParticle3r3p<Vector_3> for RelativisticParticle_3r3p {
    fn r(&self) -> &Vector_3 { &self.r }
    fn p(&self) -> &Vector_3 { &self.p }
    fn g(&self) -> &f64 { &self.g }
}*/

// just for testing, not really synchrotron...
/*struct Synchrotron_pusher { }

impl SimplePusher<RelativisticParticle_3r3p> for Synchrotron_pusher {
    fn push(&self, particle: RelativisticParticle_3r3p) -> RelativisticParticle_3r3p {
        return particle;
    }
}*/
