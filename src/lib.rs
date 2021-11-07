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

// trait for particles...
trait RelativisticParticle3r3p<V3: Vector3>: Particle {
    fn r(&self) -> &V3;
    fn p(&self) -> &V3;
    fn g(&self) -> &f64;
}

// now we define what the particles are in the memory
struct PlainRelativisticParticle3r3p { r: SimpleV3, p: SimpleV3, g: f64 }

impl Particle for PlainRelativisticParticle3r3p {}

impl RelativisticParticle3r3p<SimpleV3> for PlainRelativisticParticle3r3p {
    fn r(&self) -> &SimpleV3 { &self.r }
    fn p(&self) -> &SimpleV3 { &self.p }
    fn g(&self) -> &f64 { &self.g }
}

// just for testing, not really synchrotron...
struct SynchrotronEmissionPusher { }

impl SimplePusher<PlainRelativisticParticle3r3p> for SynchrotronEmissionPusher {
    fn push(&self, particle: &PlainRelativisticParticle3r3p) -> &PlainRelativisticParticle3r3p {
        return &particle;
    }
}
