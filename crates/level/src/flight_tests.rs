use super::*;

#[test]
fn every_shot_flies_under_every_disturbance() {
    let mut l = shipped("twin_beams");
    assert_eq!(l.flight_count(), 2);
    assert!(l.scenarios(&[]).iter().all(|s| s.field.external.is_empty()));
    l.disturbances = vec![
        Disturbance::default(),
        Disturbance {
            e: [0.5, 0.0],
            ..Disturbance::default()
        },
        Disturbance {
            waves: vec![Wave {
                amplitude: 1.0,
                direction_deg: 0.0,
                omega: 1.0,
                phase_deg: 0.0,
            }],
            ..Disturbance::default()
        },
    ];
    let scn = l.scenarios(&[]);
    assert_eq!(scn.len(), 6);
    assert_eq!(l.flight_of(4), (1, 1));
    assert_eq!(scn[4].x0, l.scenario(1, &[]).x0);
    assert!(scn[3].field.external.is_empty());
    assert_eq!(scn[4].field.external.len(), 1);
    assert!(matches!(scn[5].field.external[0], External::Wave(_)));
}
