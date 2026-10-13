import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';

const scene = new THREE.Scene();
scene.background = new THREE.Color(0x05070b);
const camera = new THREE.PerspectiveCamera(48, innerWidth / innerHeight, 0.1, 1000);
camera.position.set(0, 7, 45);
const renderer = new THREE.WebGLRenderer({ antialias: true });
renderer.setPixelRatio(Math.min(devicePixelRatio, 2));
renderer.setSize(innerWidth, innerHeight);
document.getElementById('scene').appendChild(renderer.domElement);
const controls = new OrbitControls(camera, renderer.domElement);
controls.enableDamping = true;
scene.add(new THREE.HemisphereLight(0xffffff, 0x344055, 2.0));
const light = new THREE.DirectionalLight(0xffffff, 2.4);
light.position.set(10, 15, 20); scene.add(light);
const model = new THREE.Group(); scene.add(model);
const status = document.getElementById('status');
const tooltip = document.getElementById('atom-tooltip');
const raycaster = new THREE.Raycaster();
const pointer = new THREE.Vector2(2, 2);
let pointerActive = false;
let hoveredAtom = null;
const elementNames = { C: 'Carbon', N: 'Nitrogen', O: 'Oxygen', P: 'Phosphorus', S: 'Sulfur', H: 'Hydrogen' };
const baseNames = { DA: 'Adenine', DT: 'Thymine', DG: 'Guanine', DC: 'Cytosine', A: 'Adenine', T: 'Thymine', G: 'Guanine', C: 'Cytosine' };

const color = { C: 0x939ba7, N: 0x416eff, O: 0xff454b, P: 0xffa72c, S: 0xf2d44e, H: 0xffffff };
const covalent = { C: 0.76, N: 0.71, O: 0.66, P: 1.07, S: 1.05, H: 0.31 };
const atomRadius = { C: .32, N: .33, O: .31, P: .42, S: .42, H: .20 };
const sphereGeometry = new THREE.SphereGeometry(1, 14, 10);
const bondGeometry = new THREE.CylinderGeometry(1, 1, 1, 8);
const materials = Object.fromEntries(Object.entries(color).map(([e, c]) => [e, new THREE.MeshStandardMaterial({ color: c, roughness: .48 })]));
const bondMaterial = new THREE.MeshStandardMaterial({ color: 0x9da5ad, roughness: .6 });
const strandMaterials = { A: new THREE.MeshStandardMaterial({ color: 0x13baff }), B: new THREE.MeshStandardMaterial({ color: 0xff386f }) };
let atomMeshes = [], bondMeshes = [], traceMeshes = [];
const colorMode = document.getElementById('color-mode');
const clearHighlight = document.getElementById('clear-highlight');
const selectionInfo = document.getElementById('selection-info');
const nucleotideColors = { DA: 0x34d399, DT: 0xfbbf24, DG: 0xa78bfa, DC: 0x22d3ee,
                           A: 0x34d399, T: 0xfbbf24, G: 0xa78bfa, C: 0x22d3ee };
const nucleotideMaterials = Object.fromEntries(Object.entries(nucleotideColors).map(([name, hex]) =>
    [name, new THREE.MeshStandardMaterial({color: hex, roughness: .48})]));
const dimMaterialCache = new Map();
let selectedResidue = null;
function residueKey(a) { return `${a.chain}|${a.residue}`; }
function dimMaterial(material) {
    if (!dimMaterialCache.has(material)) {
        const m = material.clone();
        m.color.multiplyScalar(.20);
        m.transparent = true;
        m.opacity = .22;
        m.depthWrite = false;
        dimMaterialCache.set(material, m);
    }
    return dimMaterialCache.get(material);
}
function refreshAtomAppearance() {
    for (const mesh of atomMeshes) {
        const a = mesh.userData.atom;
        const baseMaterial = colorMode.value === 'nucleotide'
            ? (nucleotideMaterials[a.resName] || materials[a.element])
            : materials[a.element];
        mesh.material = selectedResidue && residueKey(a) !== selectedResidue
            ? dimMaterial(baseMaterial) : baseMaterial;
    }
    for (const mesh of bondMeshes) {
        const keys = mesh.userData.residueKeys || [];
        mesh.material = selectedResidue && !keys.includes(selectedResidue)
            ? dimMaterial(bondMaterial) : bondMaterial;
    }
}
function clearSelection() {
    selectedResidue = null;
    clearHighlight.hidden = true;
    selectionInfo.textContent = 'Click any atom to highlight its nucleotide.';
    refreshAtomAppearance();
}
colorMode.addEventListener('change', refreshAtomAppearance);
clearHighlight.addEventListener('click', clearSelection);
function clearModel() { while (model.children.length) { model.remove(model.children[0]); } atomMeshes = []; bondMeshes = []; traceMeshes = []; selectedResidue = null; clearHighlight.hidden = true; selectionInfo.textContent = 'Click any atom to highlight its nucleotide.'; model.rotation.set(0, 0, 0); }
function parsePDB(pdb) {
        const atoms = [];
        for (const line of pdb.split(/\r?\n/)) {
                const record = line.slice(0, 6).trim();
                if (record !== 'ATOM' && record !== 'HETATM') continue;
                if (line.length < 54) continue;
                const alt = line[16]; if (alt && alt !== ' ' && alt !== 'A') continue;
                const element = (line.slice(76, 78).trim() || line.slice(12, 16).trim().replace(/^[0-9]+/, '')[0] || '').toUpperCase();
                if (!materials[element]) continue;
                const x = Number(line.slice(30, 38)), y = Number(line.slice(38, 46)), z = Number(line.slice(46, 54));
                if (![x, y, z].every(Number.isFinite)) continue;
                atoms.push({ name: line.slice(12, 16).trim(), element, chain: line[21], residue: line.slice(22, 27).trim(), resName: line.slice(17, 20).trim(), water: line.slice(17, 20).trim() === 'HOH', pos: new THREE.Vector3(x, y, z) });
        }
        return atoms;
}
function createBond(a, b, r = .095, mat = bondMaterial) {
        const v = new THREE.Vector3().subVectors(b, a), len = v.length();
        if (len < .01) return null;
        const mesh = new THREE.Mesh(bondGeometry, mat);
        mesh.position.copy(a).add(b).multiplyScalar(.5);
        mesh.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), v.normalize());
        mesh.scale.set(r, len, r); model.add(mesh); return mesh;
}
function residueNumber(s) { return parseInt(s, 10); }
function arePotentiallyBonded(a, b) {
        if (a.water || b.water) return false;
        if (a.chain !== b.chain) return false;
        const ra = residueNumber(a.residue), rb = residueNumber(b.residue);
        if (!Number.isFinite(ra) || !Number.isFinite(rb) || Math.abs(ra - rb) > 1) return false;
        if (ra !== rb) {
                // DNA backbone phosphodiester: O3' of previous nucleotide -> P of next nucleotide.
                return (ra + 1 === rb && a.name === "O3'" && b.name === 'P') ||
                        (rb + 1 === ra && b.name === "O3'" && a.name === 'P');
        }
        return true;
}
function renderPDB(pdb) {
        const all = parsePDB(pdb);
        const atoms = all.filter(a => !a.water);
        if (atoms.filter(a => !a.water).length < 100) throw new Error('No complete DNA atomic model found in this PDB file.');
        clearModel();
        const center = new THREE.Vector3(); const dna = atoms.filter(a => !a.water);
        dna.forEach(a => center.add(a.pos)); center.divideScalar(dna.length);
        atoms.forEach(a => a.pos.sub(center));
        // Approximate helical axis from phosphate/sugar atoms at the two ends of chain A.
        const strandA = dna.filter(a => a.chain === 'A');
        const residues = [...new Set(strandA.map(a => residueNumber(a.residue)))].sort((a, b) => a - b);
        if (residues.length >= 4) {
                const centroid = (n) => { const v = new THREE.Vector3(); const arr = strandA.filter(a => residueNumber(a.residue) === n); arr.forEach(a => v.add(a.pos)); return v.divideScalar(arr.length); };
                const axis = centroid(residues.at(-1)).sub(centroid(residues[0])).normalize();
                const q = new THREE.Quaternion().setFromUnitVectors(axis, new THREE.Vector3(1, 0, 0));
                atoms.forEach(a => a.pos.applyQuaternion(q));
        }
        for (const a of atoms) {
                const sphere = new THREE.Mesh(sphereGeometry, materials[a.element]);
                sphere.position.copy(a.pos); const r = atomRadius[a.element] || .3; sphere.scale.setScalar(a.water ? .22 : r);
                sphere.userData.atom = { element: a.element, name: a.name, chain: a.chain, residue: a.residue, resName: a.resName };
                model.add(sphere); atomMeshes.push(sphere);
        }
        // Pair search within each residue plus the known adjacent O3'-P linkage.
        // Distance-based connectivity is an approximation, not experimentally assigned bond order.
        const buckets = new Map();
        atoms.forEach((a, i) => { const k = a.chain + '|' + a.residue; if (!buckets.has(k)) buckets.set(k, []); buckets.get(k).push(i); });
        const edges = new Set();
        function tryBond(i, j) {
                const a = atoms[i], b = atoms[j]; if (!arePotentiallyBonded(a, b)) return;
                const d = a.pos.distanceTo(b.pos), limit = covalent[a.element] + covalent[b.element] + .19;
                if (d > .9 && d < limit) { const key = Math.min(i, j) + ',' + Math.max(i, j); if (!edges.has(key)) { edges.add(key); const mesh = createBond(a.pos, b.pos); if (mesh) { mesh.userData.residueKeys = [residueKey(a), residueKey(b)]; bondMeshes.push(mesh); } } }
        }
        for (const arr of buckets.values()) for (let i = 0; i < arr.length; i++)for (let j = i + 1; j < arr.length; j++)tryBond(arr[i], arr[j]);
        for (const chain of ['A', 'B']) {
                const phosphates = atoms.map((a, i) => [a, i]).filter(([a]) => a.chain === chain && a.name === 'P').sort((a, b) => residueNumber(a[0].residue) - residueNumber(b[0].residue));
                for (let i = 1; i < phosphates.length; i++) {
                        const previous = atoms.findIndex(a => a.chain === chain && residueNumber(a.residue) === residueNumber(phosphates[i][0].residue) - 1 && a.name === "O3'");
                        if (previous >= 0) tryBond(previous, phosphates[i][1]);
                }
        }
        // Backbone overlay follows sugar C1' atoms; for readability, not a covalent-bond depiction.
        for (const chain of ['A', 'B']) {
                const path = atoms.filter(a => a.chain === chain && a.name === "C1'").sort((a, b) => residueNumber(a.residue) - residueNumber(b.residue)).map(a => a.pos);
                if (path.length >= 3) { const mesh = new THREE.Mesh(new THREE.TubeGeometry(new THREE.CatmullRomCurve3(path), path.length * 14, .13, 8, false), strandMaterials[chain]); model.add(mesh); traceMeshes.push(mesh); }
        }
        traceMeshes.forEach(m => m.visible = false);
        refreshAtomAppearance();
        status.textContent = "";
        controls.target.set(0, 0, 0); controls.update();
}
async function loadRemote() { status.textContent = 'Fetching 1BNA coordinates from RCSB…'; try { const response = await fetch('https://files.rcsb.org/download/1BNA.pdb'); if (!response.ok) throw new Error('HTTP ' + response.status); const pdbText = await response.text(); renderPDB(pdbText) } catch (err) { status.textContent = 'Unable to load DNA (' + err.message + '). Check your internet connection and refresh the page.'; } }
// Use the rendered atom spheres for picking, not the bonds.
// Raycasting happens after controls/model rotation so it matches the current view.
function updateAtomHover() {
    if (!pointerActive || !atomMeshes.length) {
        tooltip.style.display = 'none';
        hoveredAtom = null;
        return;
    }
    raycaster.setFromCamera(pointer, camera);
    const hit = raycaster.intersectObjects(atomMeshes, false)[0];
    if (!hit) {
        tooltip.style.display = 'none';
        hoveredAtom = null;
        return;
    }
    const atom = hit.object.userData.atom;
    if (!atom) return;
    if (hoveredAtom !== hit.object) {
        hoveredAtom = hit.object;
        const fullName = elementNames[atom.element] || atom.element;
        const base = baseNames[atom.resName] || atom.resName;
        // Use textContent to avoid treating external PDB fields as HTML.
        tooltip.replaceChildren();
        const heading = document.createElement('strong');
        heading.textContent = `${fullName} (${atom.element})`;
        tooltip.appendChild(heading);
        for (const [label, value] of [
            ['Atom', atom.name], ['Nucleotide', base],
            ['Residue', atom.residue], ['Chain', atom.chain]
        ]) {
            const line = document.createElement('div');
            line.textContent = `${label}: ${value}`;
            tooltip.appendChild(line);
        }
    }
    tooltip.style.display = 'block';
}

renderer.domElement.addEventListener('pointermove', (event) => {
    const rect = renderer.domElement.getBoundingClientRect();
    pointer.x = ((event.clientX - rect.left) / rect.width) * 2 - 1;
    pointer.y = -((event.clientY - rect.top) / rect.height) * 2 + 1;
    pointerActive = true;
    // Keep tooltip within the viewport.
    const width = tooltip.offsetWidth || 190;
    const height = tooltip.offsetHeight || 120;
    tooltip.style.left = `${Math.max(8, Math.min(event.clientX + 16, innerWidth - width - 8))}px`;
    tooltip.style.top = `${Math.max(8, Math.min(event.clientY + 16, innerHeight - height - 8))}px`;
});
// Select an entire residue by clicking an atom; dragging to orbit does not select.
let pointerDown = null;
renderer.domElement.addEventListener('pointerdown', event => {
    if (event.button === 0) pointerDown = {x: event.clientX, y: event.clientY};
});
renderer.domElement.addEventListener('pointerup', event => {
    if (!pointerDown || event.button !== 0) return;
    const moved = Math.hypot(event.clientX - pointerDown.x, event.clientY - pointerDown.y);
    pointerDown = null;
    if (moved > 5 || !atomMeshes.length) return;
    const rect = renderer.domElement.getBoundingClientRect();
    const ndc = new THREE.Vector2(
        ((event.clientX - rect.left) / rect.width) * 2 - 1,
        -((event.clientY - rect.top) / rect.height) * 2 + 1
    );
    raycaster.setFromCamera(ndc, camera);
    const hit = raycaster.intersectObjects(atomMeshes, false)[0];
    if (!hit) { clearSelection(); return; }
    const a = hit.object.userData.atom;
    const key = residueKey(a);
    if (selectedResidue === key) { clearSelection(); return; }
    selectedResidue = key;
    clearHighlight.hidden = false;
    selectionInfo.textContent = `${baseNames[a.resName] || a.resName} · Chain ${a.chain} · Residue ${a.residue}`;
    refreshAtomAppearance();
});
renderer.domElement.addEventListener('pointercancel', () => { pointerDown = null; });
renderer.domElement.addEventListener('pointerleave', () => {
    pointerActive = false;
    tooltip.style.display = 'none';
});

let spinPaused = false;
const spinToggle = document.getElementById('spin-toggle');
spinToggle.addEventListener('click', () => {
    spinPaused = !spinPaused;
    spinToggle.textContent = spinPaused ? 'Resume Spin' : 'Pause Spin';
    spinToggle.setAttribute('aria-pressed', String(spinPaused));
});

const clock = new THREE.Clock();
function animate() {
    requestAnimationFrame(animate);
    const dt = Math.min(clock.getDelta(), .05);
    if (!spinPaused) model.rotateX(dt * 0.22);
    controls.update();
    updateAtomHover();
    renderer.render(scene, camera);
}
animate();
addEventListener('resize', () => { camera.aspect = innerWidth / innerHeight; camera.updateProjectionMatrix(); renderer.setSize(innerWidth, innerHeight) });
loadRemote();
